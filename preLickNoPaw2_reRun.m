function expe = preLickNoPaw2_reRun(Oro_data,precue, postcue)           % RV


% This function will wait until the user draw 2 rectangle (Region Of
% Interest ROI), one of which have to be centered on the rat's mouth and
% the other in a corner of the original image for background correction.

%for clarification of the function ask Roberto.

% input  = filename; is the name of the .AVI file.
% output = times; is the timestamps of CinePlex File.
%          crop; coordinates of the ROI around the Rat's mouth.
%          crop_norm; coordinates of a ROI used for background correction.
% These outputs will feed the "Licks_No_Paw_Cineplex_Filter_Noise_Removal"
% function.

%%
% precue = 4;
% postcue= 1;
pawfilter = 255;
for i = 1:length(Oro_data)
    file{i} = Oro_data(i).file(1:end-8);
end
% loop through the experimental session
for i = 1:size(file,2)
    
    K = VideoReader([file{i} '.AVI']);
    %
%     im1 = (read(K,fr(i)));
%     %
%     figure (i),imshow(im1(:,:,1));
    Times = dlmread([file{i} '.txt']);
    Crop = Oro_data(i).OroFacial.roi(1,:);
    Crop1= Oro_data(i).OroFacial.roi(2,:);
    expe(i,1).OroFacial.roi(1,:) = Crop;  %saving the coordinates of the ROI for the orofacial movement (pixel intensity changes) detection
    expe(i,1).OroFacial.roi(2,:) = Crop1;
    expe(i,1).Events.RSpout      = Oro_data(i).Events.RSpout;
    expe(i,1).Events.LSpout      = Oro_data(i).Events.LSpout;
    expe(i,1).Times              = Times;
    close all; clear crop; clear crop_norm; clear im1; clear K;
end
clear i;
%
for ii =  1:size(expe,1)
    fprintf('Processing session %0.1f\n',ii)
    
    Crop = Oro_data(ii).OroFacial.roi(1,:);
    Crop1= Oro_data(ii).OroFacial.roi(2,:);
    xright  = round(Crop(1,2))+round(Crop(1,4));
    xleft  = round(Crop(1,2));
    ytop   = round(Crop(1,1));
    ybottom = round(Crop(1,1))+round(Crop(1,3));  %
    
    %crop_norm is a crop of a blank area of the video so that background noise
    %can be subtracted from the imagecrop.
    xnright  = round(Crop1(1,2))+round(Crop1(1,4));
    xnleft   = round(Crop1(1,2));
    yntop    = round(Crop1(1,1));
    ynbottom = round(Crop1(1,1))+round(Crop1(1,3));
    %pawfilter sets the value of the pixel which is considered white enough to
    %pick up a paw in the image
    
    n = 0;
    expe(ii,1).OroFacial.framerate = 30;
    
    %duration is the number of seconds to do the analysis following the time
    expe(ii,1).OroFacial.postCS = postcue;
    
    %preduration is the time to analyze before the event
    expe(ii,1).OroFacial.preCS = precue;
    
    %edges provides the x-axis scale for creating a plot
    expe(ii,1).OroFacial.edges = -1*expe(ii,1).OroFacial.preCS: (1/expe(ii,1).OroFacial.framerate): ...
        (expe(ii,1).OroFacial.postCS - 1/expe(ii,1).OroFacial.framerate);
    
    frames    = round(expe(ii,1).OroFacial.postCS *expe(ii,1).OroFacial.framerate);
    preframes = round(expe(ii,1).OroFacial.preCS*expe(ii,1).OroFacial.framerate)  ;
    
    %This gives the total number of pixels
    pixel_number = (xright - xleft)*(ybottom - ytop);
    
    I = VideoReader([file{ii} '.AVI']);
    
    %This section of code progresses through each Event in the event structure.
    %It then cycles through each time of the individual event while
    %concurrently looping through each time in the video times
    %variable. Flag1 keeps track of where to restart the times loop after
    %an event is matched with the time variable. When the times loop surpasses a time
    %in the events variable, the first video frame following the event can be
    %pulled up at that time. m
    
    Event_Names = fieldnames(expe(ii).Events);
    for i = 1: length(Event_Names)
        
        %this checks whether the event has any trials within it
        if isempty(expe(ii,1).Events.(Event_Names{i}))
            
            %this sets Trial to an empty set and moves on to the next event
            Trial.(Event_Names{i}) = [];
            expe(ii,1).OroFacial.Trial.(Event_Names{i}) = [];
            Trial_no_out.(Event_Names{i}) = [];
            expe(ii,1).OroFacial.Trial_no_out.(Event_Names{i}) = [];
            continue
        end
        
        
        %Flag1 is a place keeper for the loop through the frames of the video.
        flag1 = 1;
        %This for loop iterates through each time in the events variable.
        for m = 1:length(expe(ii,1).Events.(Event_Names{i}))
            
            Break = 0;
            
            %This loop goes through each time in the video to find the first frame
            %following an event. Flag1 allows this loop to start where it left off
            %from the previous event. This makes the assumption that the event
            %variables are in chronological order.
            for n = flag1:length(Times)
                
                %this finds whether a video frame time is greater than the event
                %time
                if expe(ii,1).Events.(Event_Names{i})(m) <= Times(n)
                    
                    %This loop pulls up each frame from a time window around the
                    %event. Since each movement data point is a difference in pixel
                    %intensity from one frame to the next, the process starts on
                    %the 2nd frame of the prefame window. This results in one less
                    %datapoint than the frames + preframes total. This convention
                    %was used to match the method used in the earlier analysis.
                    
                    %h is the size of the filter to be used. This takes an average of the
                    %values within the filter mask (size h) for an individual pixel
                    h = ones(5,5)/25;
                    
                    %offset provides the time offset between the event and the timestamp of the
                    %first image following the event
                    Offset.(Event_Names{i})(m) = Times(n) - expe(ii,1).Events.(Event_Names{i})(m);
                    expe(ii,1).OroFacial.Offset.(Event_Names{i})(m) = Times(n) - expe(ii,1).Events.(Event_Names{i})(m);
                    
                    for k = 1: (frames + preframes)
                        
                        %This loop runs through frames + preframes times (if 1 sec pre and 2 sec
                        %post this is 90 times (30 hz sampling rate) through with the event occurring at 31)
                        %The psth is defined by im2, So the first frame following the event should
                        %be the preframes + 1 frame (or 31 in our example). This value is the first
                        %frame following the event minus the frame just before the event. For this reason, im1
                        %should be taken preframes + 1 before the event. In our example where there
                        %are 30 preframes and the event occurs just prior to frame 31, im1 should
                        %be frame zero and im2 frame 1 giving 30 im2 frames prior to the event.
                        im1 = (read(I,n - preframes - 2 + k));
                        image = imfilter(double(im1(:,:,1)), h);
                        %image = imrotate(image, -90);
                        
                        im2 = (read(I,(n - preframes - 1 + k)));
                        image2 = imfilter(double(im2(:,:,1)), h);
                        %image2 = imrotate(image2, -90);
                        clear im1 im2
                        
                        %inorm and val95 find aspects of the noise which are used to filter the
                        %image. The norm_Crop is used to assess difference in pixels frame to
                        %frame in an area where no intensity difference should be occurring. A
                        %histogram of the differences frame to frame (not absolute) should ideally
                        %be centered around 0 with a normal distribution of the pixel difference frame to frame.
                        %inorm is the offset of the mean of the distribution from zero. This is
                        %subtracted from the difference between the two images. Once the offset is
                        %removed, a %95 confidence interval of the noise is determined. This is set
                        %to val95. Since any changes of mouth movement that are within the frame to
                        %frame subtracted pixel noise can not be distinguished from the noise, all
                        %changes in pixel intensity below the 95% confidence value (val95) are
                        %removed and the difference matrix is reduced by this amount with any
                        %negative values set to zero.
                        %This code finds the difference in pixel activity in the norm_crop area.
                        %The difference in the norm_crop area is attributed to temporal visual noise and is
                        %subtracted from the absolute difference in the output.
                        ioffset = mean(mean(image2(xnleft:xnright, yntop: ynbottom) - image(xnleft:xnright, yntop: ynbottom)));
                        %this creates a matrix of the difference frame to frame of the Crop_norm
                        %region with the offset subtracted so as to make the mean of the
                        %distribution zero. The absolute value is taken so that only one side of
                        %the confidence interval needs to be taken (the noise is assumed to be
                        %gaussian distributed.
                        norm = abs(image2(xnleft:xnright, yntop: ynbottom) - image(xnleft:xnright, yntop: ynbottom) - ioffset);
                        
                        %this finds the difference value which includes 95% of the noise
                        norm_sorted = sort(norm);
                        
                        %To make sure there are not massive fluctuations in val95, it is averaged
                        %over the previous 5 trials
                        if k <= 5
                            vals(k) = norm_sorted(ceil(0.95*length(norm_sorted)));
                            val95 = mean(vals);
                        else
                            %circshift permutes the values in vals so that the oldest value moves
                            %to the five position and is overwritten
                            vals = circshift(vals, [1,-1]);
                            vals(5) = norm_sorted(ceil(0.95*length(norm_sorted)));
                            val95 = mean(vals);
                        end
                        
                        %image3 stores the values of the absolute difference
                        image3 = zeros((xright - xleft), (ybottom - ytop));
                        
                        %the following two loops iterate through each pixel in the image
                        %matrix, but only for the cropped area.
                        for x = 1:(xright - xleft)
                            
                            %The following if statement allows this for loop to break if a
                            %break is called in the next for loop. Break is the variable
                            %that keeps track of this
                            
                            if Break == 1
                                
                                Break = 0;
                                break
                            end
                            
                            for y = 1:(ybottom - ytop)
                                
                                %The following statements check whether the pixel and
                                %surrounding pixels are white, if so, the frame is dropped
                                %from the analysis, as a paw is probably in the frame.
                                %Pawfilter is the variable that must be set as a threshold
                                %for paw detection
                                
                                if (x > 1 && y > 1) && (((image(xleft + x - 1,ytop + y - 1,1) >= pawfilter && ...
                                        image(xleft + x - 2,ytop + y - 1,1) >= pawfilter) || ...
                                        (image(xleft + x - 1,ytop + y - 1,1) && ...
                                        image(xleft + x - 1,ytop + y - 2,1) >= pawfilter)) || ...
                                        ((image2(xleft + x - 1,ytop + y - 1,1) >= pawfilter && ...
                                        image2(xleft + x - 2,ytop + y - 1,1) >= pawfilter) ...
                                        || (image2(xleft + x - 1,ytop + y - 1,1) >= pawfilter && ...
                                        image2(xleft + x - 1,ytop + y - 2,1) >= pawfilter)))
                                    
                                    image3 = NaN((xright - xleft), (ybottom - ytop));
                                    
                                    Break = 1;
                                    break
                                end
                                
                                %Each change in pixel intensity is stored in image3, inorm is subtracted
                                %off here for each pixel change as it is the estimate for the temporal noise of the
                                %video.
                                image3(x,y) = abs(image2(xleft + x - 1,ytop + y - 1,1) - image(xleft + x - 1,ytop + y - 1,1) - ioffset);
                                
                                
                                
                            end
                        end
                        
                        %the following two lines removes all values within the noise band by first
                        %subtracting off the value marking the 95% CDF of the noise and setting all
                        %negative values to zero.
                        image3 = image3 - val95;
                        
                        image3(image3 < 0) = 0;
                        
                        %The total change in pixel intensity for a single frame is stored
                        %in the output matrix trial, which contains trials in rows and
                        %times of the frames in columns. This is normalized by the number of
                        %pixels so the total is an average per pixel
                        Trial.(Event_Names{i})(m,k) = sum(sum(image3))/pixel_number;
                        expe(ii,1).OroFacial.Trial.(Event_Names{i})(m,k) = sum(sum(image3))/pixel_number;
                        
                        clear image2 image image3;
                    end
                    
                    %Flag1 needs to be set so that the loop doesn't start all the
                    %way back at the beginning of the times variable.
                    flag1 = n;
                    
                    %This ends the n for loop since the event was found
                    break
                    
                end
            end
            
        end
        
        %The following lines remove any values in trials lying three standard
        %deviations from the mean by replacing the value with NaN and output this
        %matrix as trial_no_out
        trial_list = expe(ii,1).OroFacial.Trial.(Event_Names{i})(:);
        meantrial = nanmean(trial_list);
        stdtrial = nanstd(trial_list);
        thresh_high = meantrial + 1*stdtrial;
        [b1, b2] = find(expe(ii,1).OroFacial.Trial.(Event_Names{i}) >= thresh_high);
        expe(ii,1).OroFacial.Trial_no_out.(Event_Names{i}) = expe(ii,1).OroFacial.Trial.(Event_Names{i});
        
        for j = 1 : length(b1)
            expe(ii,1).OroFacial.Trial_no_out.(Event_Names{i})(b1(j), b2(j)) = NaN;
        end
        
        clear trial_list meantrial stdtrial thresh_high b1 b2
        
    end
    expe(ii,1).auROC.RSpout = mean(expe(ii,1).OroFacial.Trial.RSpout)';
    expe(ii,1).auROC.LSpout = mean(expe(ii,1).OroFacial.Trial.LSpout)';

end

