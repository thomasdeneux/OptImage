function x = fn4D_polyunion(varargin)

x = zeros(2,1);
for i=1:length(varargin)
    if ~isempty(varargin{i});
        x = [x NaN(2,1) varargin{i}];
    end
end
