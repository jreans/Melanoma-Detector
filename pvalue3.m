function [M2,ptable] = pvalue3(X,y)
ptable = zeros(1,size(X,2));
melpos = find(y == 1);
mel = X(melpos,:);
nevpos = find(y == 0);
nev = X(nevpos,:)
for i_im = 1:size(X,2)
    [~,p1] = ttest2(mel(:,i_im), nev(:,i_im));
    ptable(1,i_im) = p1;
end
M2 = find(ptable <= 0.05);
end