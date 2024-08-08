%%PCA analysis of ground-truth data from analysis of simulated actin networks
simulated_data = [reshape(LFB,450,1),reshape(all_occup_per_site,450,1),...
    reshape(order_parameter_in_network,450,1),reshape(mean_angles_in_network,450,1),...
    reshape(angular_variation_in_network,450,1),...
    reshape(all_median_per_site,450,1)]; %LFB|Occupancy|Order Parameter|Mean Angle|Angular Variation|Distance

end_of_array = size(simulated_data,2)+1;
for i = 1:9
    xlink = (i-1)*100;
    start = (i-1)*50+1;
    for j = start:start+49
        simulated_data(j,end_of_array) = xlink;
    end
end

simulated_data = array2table(simulated_data,...
    'VariableNames',{'LFB','Occupancy','Order Parameter','Mean Angle','Angular Variation','Distance','Nc'});

numericVars_simulated = varfun(@isnumeric,simulated_data,"output","uniform");
numericPredictorTable_simulated = simulated_data(:, numericVars_simulated);
numericPredictorMatrix_simulated = table2array(varfun(@double, numericPredictorTable_simulated)); %total_synthetic_data;%
[coeffs_simulated, transformedData_simulated, ~, ~, explained_simulated] = pca(zscore(numericPredictorMatrix_simulated)); %zscore();
enoughExplained = cumsum(explained_simulated)/sum(explained_simulated) >= 95/100;
numberOfComponentsToKeep_simulated = find(enoughExplained, 1);
disp("The number of components needed to explain at least " + 95 + "% of the variance is " + num2str(numberOfComponentsToKeep_simulated))

figure
h = biplot(coeffs_simulated(:,[1 2]),'scores',transformedData_simulated(:,[1 2]));
hID = get(h, 'tag');
hPt = h(strcmp(hID,'obsmarker'));
hLine = h(strcmp(hID,'varline'));
hVar = h(strcmp(hID,'varmarker')); 
grp = table2array(simulated_data.Nc);
grpID = unique(grp);

for i = 1:length(grpID)
set(hPt(grp==grpID(i)), 'Color', color_cells{i,1},'MarkerSize',100);%, 'DisplayName', sprintf('Cluster %d', grpID(i)))
end
for i = 1:length(hLine)
set(hLine(i),'LineStyle','none','Marker','none');
set(hVar(i),'Marker','none');
end
hold on
g = biplot(coeffs_simulated(:,[1 2]),"varLabel",simulated_data.Properties.VariableNames);
gID = get(g, 'tag');
gtext = g(strcmp(gID,'varlabel'));
gLine = g(strcmp(gID,'varline'));
for i = 1:length(gtext)
set(gtext(i),'FontSize',20,'FontWeight','Bold');
end
for i = 1:length(gLine)
set(gLine(i),'Linewidth',4,'Marker','none');
end
set(gcf,'color','w');
set(gca,'box','on','LineWidth',5,'color','w','Fontsize',80,'FontName','Helvetica Neue');
pbaspect([1 1 1]);


%%PCA analysis of morphometric data from image analysis of pseudo-fluorescence images
Morpho_from_images = [];
for i = 1:9
    num = i-1;
    name  = strcat('Morphometrycorrected2v17bS',num2str(num)); %change filename to match output from ImageJ macro
    bundling_level = eval(name);
    bundling_level.cvDist = bundling_level.CoefficientOfVariation.*bundling_level.MedianDistanceFF;
    count = 1;
    p = table2array(bundling_level(1:50,[33,7,20,17,18,28])); %Bundle Parameter|Occupancy|Order Parameter|Mean Angle|Angular Variation|Distance
    Morpho_from_images = cat(1,Morpho_from_images,p);
end

end_of_array = size(Morpho_from_images,2)+1;
for i = 1:9
    xlink = (i-1)*100;
    start = (i-1)*50+1;
    for j = start:start+49
        Morpho_from_images(j,end_of_array) = xlink;
    end
end

Morpho_from_images=array2table(Morpho_from_images,...
    'VariableNames',{'Bundle Parameter','Occupancy','Order Parameter','Mean Angle','Angular Variation','Distance','Nc'});

numericVars_images = varfun(@isnumeric,Morpho_from_images,"output","uniform");
numericPredictorTable_images = Morpho_from_images(:, numericVars_images);
numericPredictorMatrix_images = table2array(varfun(@double, numericPredictorTable_images)); %total_synthetic_data;%
[coeffs_images, transformedData_images, ~, ~, explained_images] = pca(zscore(numericPredictorMatrix_images)); %feature normalization followed by pca dimensionality reduction
enoughExplained = cumsum(explained_images)/sum(explained_images) >= 95/100;
numberOfComponentsToKeep_images = find(enoughExplained, 1);
disp("The number of components needed to explain at least " + 95 + "% of the variance is " + num2str(numberOfComponentsToKeep_images))

figure
h = biplot(coeffs_images(:,[1 2]),'scores',transformedData_images(:,[1 2]));
hID = get(h, 'tag');
hPt = h(strcmp(hID,'obsmarker'));
hLine = h(strcmp(hID,'varline'));
hVar = h(strcmp(hID,'varmarker')); 

grp = table2array(pca_table3(:,8));
grpID = unique(grp);

for i = 1:length(grpID)
set(hPt(grp==grpID(i)), 'Color', color_cells{i,1},'MarkerSize',100);
end
for i = 1:length(hLine)
set(hLine(i),'LineStyle','none','Marker','none');
set(hVar(i),'Marker','none');
end
hold on
g = biplot(coeffs_simulated(:,[1 2]),"varLabel",Morpho_from_images.Properties.VariableNames);
gID = get(g, 'tag');
gtext = g(strcmp(gID,'varlabel'));
gLine = g(strcmp(gID,'varline'));
for i = 1:length(gtext)
set(gtext(i),'FontSize',20,'FontWeight','Bold');
end
for i = 1:length(gLine)
set(gLine(i),'Linewidth',4,'Marker','none');
end
set(gcf,'color','w');
set(gca,'box','on','LineWidth',5,'color','w','Fontsize',80,'FontName','Helvetica Neue');
pbaspect([1 1 1]);
