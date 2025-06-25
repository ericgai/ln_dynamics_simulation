fid = fopen('density_counts.txt');
lines = textscan(fid,'%s','delimiter','\n');
fclose(fid);
lines = lines{1};
ova_arr = str2double ( split(lines{1},',') );
klh_arr = str2double ( split(lines{2},',') );
sim_arr = str2double ( split(lines{3},',') );

fig = figure;

ax1=subplot(3,1,1);
h1 = histogram(ova_arr,norm='probability');
xticks([]);
title('weak immune response (experiment)')
%xlim([-1 24]);

ax2=subplot(3,1,2);
h2 = histogram(klh_arr,norm='probability');
xticks([]);
title('strong immune response (experiment)')

%xlim([-1 24]);


ax3=subplot(3,1,3);
h3 = histogram(sim_arr,norm='probability');
title('strong immune response (simulation)')

linkaxes([ax1,ax2,ax3],'xy');

han=axes(fig,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'probability density');
xlabel(han,'activation density around DCs');

set(1,'renderer','painters');
