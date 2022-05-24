function [peaks_new_val,peaks_new_pos]=RealPeaks(peak_pos,sol,t,ton,Threshold)
    
    peaks = peak_pos(t(peak_pos)>ton);
    peaks_new_pos = []; n = 1;
    peaks_new_val = [];

    if length(peaks)==1
        if and((sol(peaks)/min(sol(1:peaks)))>=Threshold,(sol(peaks)/min(sol(peaks:end)))>=Threshold)
            peaks_new_pos(n) = peaks;
            peaks_new_val(n) = sol(peaks);
            n=n+1;
        end
    else
        
        for j = 1:length(peaks);
            if j == 1;
                if and((sol(peaks(j))/min(sol(1:peaks(j))))>=Threshold,(sol(peaks(j))/min(sol(peaks(j):peaks(j+1))))>=Threshold)
                    peaks_new_pos(n) = peaks(j);
                    peaks_new_val(n) = sol(peaks(j));
                    n=n+1;
                end
            elseif j == length(peaks) % Peak has to be significantly higher than surroundings and return below a certain value.
                if (sol(peaks(j))/min(sol(peaks(j):end)))>=Threshold
                    peaks_new_pos(n) = peaks(j);
                    peaks_new_val(n) = sol(peaks(j));
                    n=n+1;
                end
            else
                if and((sol(peaks(j))/min(sol(peaks(j):peaks(j+1))))>=Threshold,(sol(peaks(j))/min(sol(peaks(j-1):peaks(j))))>=Threshold)
                    peaks_new_pos(n) = peaks(j);
                    peaks_new_val(n) = sol(peaks(j));
                    n=n+1;
                end
            end
        end
    end
end