%*************************************************************************%
% Project: Function to Split Vector into multiple categories/components
% Author: James McCanna
% Date Started: 7/7/14
% Date Last Modified: 4/13/16
% Filename: splitVectorCat.m
%*************************************************************************%

function [output] = splitVectorCat(input,index)

categories = unique(index);
colCount = zeros(1,length(categories));
colCount(colCount == 0) = 1;
output = NaN(length(input),length(colCount));

for counter = 1:1:length(input)
    for col = 1:1:length(colCount)
        %if col > 1
        if isnumeric(categories(col))
            if index(counter) == categories(col)
                output(colCount(col),col) = input(counter);
                colCount(col) = colCount(col)+1;
            end
        else
            if strcmp(index(counter),categories(col))
                output(colCount(col),col) = input(counter);
                colCount(col) = colCount(col)+1;
            end
        end
        %elseif index(counter) < upperLimits(col)
        %output(colCount(col),col) = input(counter);
        %colCount(col) = colCount(col)+1;
        %end
    end
    %     if index(counter) < upperLimits(1)
    %         output(colCount(1),1) = input(counter);
    %         colCount(1) = colCount(1)+1;
    %     elseif length(upperLimits) > 2 && index(counter) < upperLimits(2)
    %         output(colCount(2),2) = input(counter);
    %         colCount(2) = colCount(2)+1;
    %     elseif length(upperLimits) > 3 && index(counter) < upperLimits(3)
    %         output(colCount(3),3) = input(counter);
    %         colCount(3) = colCount(3)+1;
    %     elseif length(upperLimits) > 4 && index(counter) < upperLimits(4)
    %         output(colCount(4),4) = input(counter);
    %         colCount(4) = colCount(4)+1;
    %     elseif length(upperLimits) > 5 && index(counter) < upperLimits(5)
    %         output(colCount(5),5) = input(counter);
    %         colCount(5) = colCount(5)+1;
    %     else
    %         output(colCount(length(colCount)),length(colCount)) = input(counter);
    %         colCount(length(colCount)) = colCount(length(colCount))+1;
    %     end
end

% Remove NaN values from output
%for i=1:1:size(output,2)
%    output(isnan(output(:,i)),:)=[];
%end

end