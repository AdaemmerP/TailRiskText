function Output = check_vintage_series(search_var, T, dates)
present = true;
i = 1;
while present && i <= T 
    disp([num2str(i*100/T) '% of vintages checked'])
    [~, header2] = xlsread([dates(i,:), '.csv']); %#ok<*XLSRD> 
    var_names = header2(1,2:end);
    select = contains(var_names, search_var);
    if sum(select) > 0
        break
    end
    i = i + 1;
end
Output = i;
