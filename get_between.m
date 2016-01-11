function out = get_between(str, s1, s2);

dum = regexp(yy, str , 'split')
count = 1;
for i = 2:length(dum)
    dum2 = regexp(dum{i}, '</a>', 'split');
    data(count).date = dum2{1};
    count = count + 1;
end