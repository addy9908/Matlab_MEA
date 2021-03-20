%% Author: Zengyou Ye
% credit: adapt from online source
% goal: change the column number to column letter, 1st row in excel
% example: number 2 to B1
function colName = num2col(colNum)
    % create excel number to column
    z = 'A':'Z';
    colLetter = arrayfun(@(x)z(rem(floor(x*26.^(1-floor(log(x)/log(26)+1):0)),26)),colNum,'un',0);
    colName = sprintf('%s1', colLetter{1});
end