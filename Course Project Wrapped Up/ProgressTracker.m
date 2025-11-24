classdef ProgressTracker < handle
    properties
        Count = 0;
        Total = 1;
        StartTime
    end

    methods
        function obj = ProgressTracker(total)
            obj.Total = total;
            obj.StartTime = tic;
        end

        function update(obj)
            obj.Count = obj.Count + 1;
            if mod(obj.Count, max(1, floor(obj.Total / 1000))) == 0 || obj.Count == obj.Total
                elapsed = toc(obj.StartTime);
                remaining = (elapsed / obj.Count) * (obj.Total - obj.Count);
                fprintf('Progress: %5.1f%% (%d / %d), Elapsed: %.1fs, Remaining: %.1fs\n', ...
                        100 * obj.Count / obj.Total, obj.Count, obj.Total, elapsed, remaining);
            end
        end
    end
end
