%% image analysis
% with some linux based commands
clc
clear
close all

% use bioformats
addpath(genpath('~/Documents/MATLAB/shared/bfmatlab'))

% create output folder
!mkdir -p output

%% image analysis

% coffee ring radius
r_coffee = 1.5316e+03;
% radius range where the pattern will be evaluated
R_ratio = [.875,.95];
% radius step size in calculations
dr = 1; % / pixel

% file names
names = ['A':'H','W'];
% file numbers
numbers = 1:6;

% memory allocation for the given properties
s_red_acf = zeros(numel(numbers),numel(names));
s_green_acf = zeros(numel(numbers),numel(names));
redratio = zeros(numel(numbers),numel(names));
width = zeros(numel(numbers),numel(names));

hwait = waitbar(0,'Please wait...');  % Initialize waitbar
tic
for name_act = 1:numel(names)
    name = names(name_act);

    % name of the folder for the images
    [~,path_names] = system(['find ~+/2023-08-02_ABCDEFGH_WT_Day5_pEP28_pEp17/' name '*/ -name ''*.czi''']);

    path_names = splitlines(strtrim(path_names));

    if isempty(path_names)
        error('Unnown image.')
    end

    serie = 1;

    for ssz = numbers

        path_name = path_names(ssz);
        % the data file with the processed image and the fitted ellipse
        mat_name = strrep(path_name{1},'czi','mat');
        % variables in the file
        if exist(mat_name,'file')
            variables = who('-file', mat_name);
        else
            variables = [];
        end


        % caluclate matrix from the max values
        % scale: um/pixel
        if ismember('Cmax', variables)
            load(mat_name,'Cmax','scale')
        else
            disp('Process image...')
            tic
            [Cmax,scale] = draw_max(path_name,serie);
            toc
            save(mat_name,'Cmax','scale')
        end

        A = zeros(size(Cmax,1),size(Cmax,2),3);

        thresh = multithresh(Cmax,1);
        A(:,:,1) = imquantize(Cmax(:,:,1),thresh);

        thresh = multithresh(Cmax,2);
        A(:,:,2) = imquantize(Cmax(:,:,2),thresh);

        % we did not save the ellipse before
        if ~isempty(variables) && (ismember('stats', variables) || ismember('Rmax', variables))
            load(mat_name,'stats','Rmax')
        else
            tic
            disp('Find ellipse...')
            [stats,Rmax] = find_ellipse(A);
            toc
            % visellipse(stats,r_coffee/scale,'b')
            save(mat_name,'stats','Rmax','-append')
        end

        R1 = Rmax*R_ratio(1);
        R2 = Rmax*R_ratio(2);

        % Calculation of the parameters
        [sep_red_acf,sep_green_acf,ratio] = calc_polar(Cmax,dr,R1,R2,stats,scale);

        % store the data
        s_red_acf(ssz,name_act) = mean(sep_red_acf);
        s_green_acf(ssz,name_act) = mean(sep_green_acf);
        redratio(ssz,name_act) = mean(ratio);
        width(ssz,name_act) = Rmax*scale-r_coffee;

    end

    % Update the waitbar
    if mod(name_act,max(round(numel(names)/100),1))==0
        waitbar(name_act / numel(names), hwait, ['Progress: ', int2str(name_act/numel(names)*100), '%']);
    end
end
toc
close(hwait);  % Close the waitbar

%% bisquare weights

% memory allocation for the different properties
for prop = ["green","red","width","redratio"]
    weighted_mean.(prop) = zeros(1,numel(names.'));
    confint_calced.(prop) = zeros(1,numel(names.'));
    weights.(prop) = zeros(size(s_red_acf));
end

% calculation
for i = 1:numel(names.')
    [weighted_mean.green(1,i),confint_calced.green(1,i),weights.green(:,i)] = stat_calc(s_green_acf(:,i));
    [weighted_mean.red(1,i),confint_calced.red(1,i),weights.red(:,i)] = stat_calc(s_red_acf(:,i));
    [weighted_mean.width(1,i),confint_calced.width(1,i),weights.width(:,i)] = stat_calc(width(:,i));
    [weighted_mean.redratio(1,i),confint_calced.redratio(1,i),weights.redratio(:,i)] = stat_calc(redratio(:,i));
end

%% plot

figure
s_red2 = s_red_acf;
s_red2(weights.red==0) = nan;
boxplot(s_red2,'Labels',names.')
ylabel('s_{red,acf}')
exportgraphics(gcf,'output/red_acf.pdf')

figure
s_green2 = s_green_acf;
s_green2(weights.red==0) = nan;
boxplot(s_green2,'Labels',names.')
ylabel('s_{green,acf}')
exportgraphics(gcf,'output/green_acf.pdf')

figure
tmp = redratio;
tmp(weights.redratio==0) = nan;
boxplot(tmp,'Labels',names.')
ylabel('x_{red}')
exportgraphics(gcf,'output/redratio.pdf')

figure
tmp = width;
tmp(weights.width==0) = nan;
boxplot(tmp,'Labels',names.')
ylabel('width')
exportgraphics(gcf,'output/width.pdf')

target = weighted_mean.green./weighted_mean.red;
disp(target(1:2:end-1)./target(2:2:end))

target = mean(s_green2,"omitnan")./mean(s_red2,"omitnan");
disp(target(1:2:end-1)./target(2:2:end))


%% export data

readme = cell(1,5);
readme(1,1) = {'name'};
readme(1,2) = {'red ratio'};
readme(1,3) = {'width / um'};
readme(1,4) = {'s(red) / um'};
readme(1,5) = {'s(green) / um'};
!rm output/data.xls
writecell(readme,'output/data.xls','Sheet','data')
writecell(readme,'output/data.xls','Sheet','outliners')
writecell(readme,'output/data.xls','Sheet','weighted_mean')

% first line where we can write
act_line = 2;

% data number
n = numel(names)*numel(numbers);

% names for the lines
name_array = char(ones(numel(numbers),1)*names);
name_array = [name_array(:),int2str(reshape(numbers.'*ones(1,numel(names)),numel(numbers)*numel(names),1))];

% write to xls file
writematrix(name_array,'output/data.xls','Sheet','data','Range',['A' int2str(act_line) ':A' int2str(act_line+n-1)])
writematrix(redratio(:),'output/data.xls','Sheet','data','Range',['B' int2str(act_line) ':B' int2str(act_line+n-1)])
writematrix(width(:),'output/data.xls','Sheet','data','Range',['C' int2str(act_line) ':C' int2str(act_line+n-1)])
writematrix(s_red_acf(:),'output/data.xls','Sheet','data','Range',['D' int2str(act_line) ':D' int2str(act_line+n-1)])
writematrix(s_green_acf(:),'output/data.xls','Sheet','data','Range',['E' int2str(act_line) ':E' int2str(act_line+n-1)])


% outliners
writematrix(name_array,'output/data.xls','Sheet','outliners','Range',['A' int2str(act_line) ':A' int2str(act_line+n-1)])
writematrix(weights.redratio(:),'output/data.xls','Sheet','outliners','Range',['B' int2str(act_line) ':B' int2str(act_line+n-1)])
writematrix(weights.width(:),'output/data.xls','Sheet','outliners','Range',['C' int2str(act_line) ':C' int2str(act_line+n-1)])
writematrix(weights.red(:),'output/data.xls','Sheet','outliners','Range',['D' int2str(act_line) ':D' int2str(act_line+n-1)])
writematrix(weights.green(:),'output/data.xls','Sheet','outliners','Range',['E' int2str(act_line) ':E' int2str(act_line+n-1)])

% weighted mean
writematrix(names.','output/data.xls','Sheet','weighted_mean','Range',['A' int2str(act_line) ':A' int2str(act_line+numel(names)-1)])
writematrix(weighted_mean.redratio(:),'output/data.xls','Sheet','weighted_mean','Range',['B' int2str(act_line) ':B' int2str(act_line+numel(names)-1)])
writematrix(weighted_mean.width(:),'output/data.xls','Sheet','weighted_mean','Range',['C' int2str(act_line) ':C' int2str(act_line+numel(names)-1)])
writematrix(weighted_mean.red(:),'output/data.xls','Sheet','weighted_mean','Range',['D' int2str(act_line) ':D' int2str(act_line+numel(names)-1)])
writematrix(weighted_mean.green(:),'output/data.xls','Sheet','weighted_mean','Range',['E' int2str(act_line) ':E' int2str(act_line+numel(names)-1)])