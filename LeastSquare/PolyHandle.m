function [f,fstr] = PolyHandle(a)

deg = length(a) - 1;

if a(end) ~= 0
    fstr = ['(',num2str(a(end)),')*x^',num2str(deg)];
else
    fstr = '';
end

for ii = 1:(deg-1)
    if a(end-ii) ~= 0
        if ~isempty(fstr)
            fstr = strcat(fstr,' + (',num2str(a(end-ii)),')*x^',num2str(deg-ii));
        else
            fstr = strcat('(',num2str(a(end-ii)),')*x^',num2str(deg-ii));
        end
    end
end

fstr = strcat(fstr,' + (',num2str(a(1)),')');
fstr = replace(fstr,'x^1','x');

f = eval(['@(x) ', replace(fstr,'^','.^')]);