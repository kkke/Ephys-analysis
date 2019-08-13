function Sum = tsTasteError(Sum)
for i = 1:length(Sum)
    M = Sum(i).event.tsRCorr.TasteID.M;
    O = Sum(i).event.tsRCorr.TasteID.O;
    TasteDel = Sum(i).event.tsRCorr.TasteDel;
    TasteID_r = unique(TasteDel(:,2));
    r1 = TasteDel(find(TasteDel(:,2)==TasteID_r(1)),:);
    r2 = TasteDel(find(TasteDel(:,2)==TasteID_r(2)),:);
    if size(r1,1) == size(M,1) && size(r2,1) == size(O,1)&& size(r2,1) == size(M,1) && size(r1,1) == size(O,1);
        ts_dif1 = max(M-r1(:,1));
        ts_dif2 = max(O-r2(:,1));
        if ts_dif1<2 && ts_dif1>1 && ts_dif2<2 && ts_dif2>1
            TasteID_rT{1} = 'M';
            TasteID_rT{2} = 'O';
        else
            ts_dif1 = max(M-r2(:,1));
            ts_dif2 = max(O-r1(:,1));
            if ts_dif1<2 && ts_dif1>1 && ts_dif2<2 && ts_dif2>1
                TasteID_rT{1} = 'O';
                TasteID_rT{2} = 'M';
            else
                fprintf('%f',i)
                error('Something is wrong')
            end
        end
    elseif size(r1,1) == size(M,1) && size(r2,1) == size(O,1);
        ts_dif1 = max(M-r1(:,1));
        ts_dif2 = max(O-r2(:,1));
        if ts_dif1<2 && ts_dif1>1 && ts_dif2<2 && ts_dif2>1
            TasteID_rT{1} = 'M';
            TasteID_rT{2} = 'O';
        else
            fprintf('%f',i)
            error('Something is wrong')
        end
    elseif size(r2,1) == size(M,1) && size(r1,1) == size(O,1);
        ts_dif1 = max(M-r2(:,1));
        ts_dif2 = max(O-r1(:,1));
        if ts_dif1<2 && ts_dif1>1 && ts_dif2<2 && ts_dif2>1
            TasteID_rT{1} = 'O';
            TasteID_rT{2} = 'M';
        else
            fprintf('%f',i)
            error('Something is wrong')
        end
        
    end
    if size(Sum(i).event.tsRErr.FLickTaste,1) == size(Sum(i).event.tsRErr.TasteDel,1)
        FLickTaste = [Sum(i).event.tsRErr.FLickTaste,Sum(i).event.tsRErr.TasteDel(:,2)];
        Sum(i).event.tsRErr.TasteID.(TasteID_rT{1}) = FLickTaste(find(FLickTaste(:,2)==TasteID_r(1)),1);
        Sum(i).event.tsRErr.TasteID.(TasteID_rT{2}) = FLickTaste(find(FLickTaste(:,2)==TasteID_r(2)),1);
    else
        fprintf('%f',i)
        error('Something is wrong')
    end 
end
%%
for i = 1:length(Sum)
    S = Sum(i).event.tsLCorr.TasteID.S;
    Q = Sum(i).event.tsLCorr.TasteID.Q;
    TasteDel = Sum(i).event.tsLCorr.TasteDel;
    TasteID_r = unique(TasteDel(:,2));
    r1 = TasteDel(find(TasteDel(:,2)==TasteID_r(1)),:);
    r2 = TasteDel(find(TasteDel(:,2)==TasteID_r(2)),:);
    if size(r1,1) == size(S,1) && size(r2,1) == size(Q,1)&& size(r2,1) == size(S,1) && size(r1,1) == size(Q,1);
        ts_dif1 = max(S-r1(:,1));
        ts_dif2 = max(Q-r2(:,1));
        if ts_dif1<2 && ts_dif1>1 && ts_dif2<2 && ts_dif2>1
            TasteID_rT{1} = 'S';
            TasteID_rT{2} = 'Q';
        else
            ts_dif1 = max(S-r2(:,1));
            ts_dif2 = max(Q-r1(:,1));
            if ts_dif1<2 && ts_dif1>1 && ts_dif2<2 && ts_dif2>1
                TasteID_rT{1} = 'Q';
                TasteID_rT{2} = 'S';
            else
                fprintf('%f',i)
                error('Something is wrong')
            end
        end
    elseif size(r1,1) == size(S,1) && size(r2,1) == size(Q,1);
        ts_dif1 = max(S-r1(:,1));
        ts_dif2 = max(Q-r2(:,1));
        if ts_dif1<2 && ts_dif1>1 && ts_dif2<2 && ts_dif2>1
            TasteID_rT{1} = 'S';
            TasteID_rT{2} = 'Q';
        else
            fprintf('%f',i)
            error('Something is wrong')
        end
    elseif size(r2,1) == size(S,1) && size(r1,1) == size(Q,1);
        ts_dif1 = max(S-r2(:,1));
        ts_dif2 = max(Q-r1(:,1));
        if ts_dif1<2 && ts_dif1>1 && ts_dif2<2 && ts_dif2>1
            TasteID_rT{1} = 'Q';
            TasteID_rT{2} = 'S';
        else
            fprintf('%f',i)
            error('Something is wrong')
        end
        
    end
    if size(Sum(i).event.tsLErr.FLickTaste,1) == size(Sum(i).event.tsLErr.TasteDel,1)
        FLickTaste = [Sum(i).event.tsLErr.FLickTaste,Sum(i).event.tsLErr.TasteDel(:,2)];
        Sum(i).event.tsLErr.TasteID.(TasteID_rT{1}) = FLickTaste(find(FLickTaste(:,2)==TasteID_r(1)),1);
        Sum(i).event.tsLErr.TasteID.(TasteID_rT{2}) = FLickTaste(find(FLickTaste(:,2)==TasteID_r(2)),1);
    else
        fprintf('%f',i)
        error('Something is wrong')
    end 
end
end