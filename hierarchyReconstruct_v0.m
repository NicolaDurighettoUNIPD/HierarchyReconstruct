clear;
close all;
clc;
% load('viterbo_all.mat');
load('valfredda.mat');

%% restructure data and create graph

% Definition of observation dates. This function can be used to get dates from available data, or it can be done manually.
dates = getDates(network);
% dates = datetime(2020, 01, 01):datetime(2020, 01, 12);

% Node names, used almost only for plots
nodeId = network.Name;

% Matrix of observed states: one row for each day, one column for each node
statusObs = getStatusMatrix(network, dates);

% Create the first graph, O. Here, node A is connected to B only if there is an observation of wet A and dry B (i.e. a clue that A is more persistent than B). The weight associated to edge A->B is the number of such observations. There can also be and edge B->A, if there are such observations in the dataset. 
observationMatrix = getObservationMatrix(statusObs);

% Create graph H, describing the hierarcy. Here, loops are removed.
hierarchyMatrix = getHierarchyMatrix(observationMatrix);
hierarchyMatrix = removeLoops(hierarchyMatrix);
hierarchyGraph = digraph(hierarchyMatrix, network);

%% Reorder graph topologically
order = toposort(hierarchyGraph);

network = network(order, :);
nodeId = nodeId(order);
statusObs = statusObs(:, order);
observationMatrix = observationMatrix(order, order);
hierarchyMatrix = hierarchyMatrix(order, order);

%% Check hierarchy accuracy
% Option 1: calculate the fraction of pairwise observations in accordance with the hierarcy. Since we reordered the matrix, these are all in the lower triangle.

hierarchyAccuracy = sum(tril(observationMatrix), 'all') / sum(observationMatrix, 'all');

% Optrion 2: use the hierarchy to reconstruct the status of the nodes, and calculate the corresponding accuracy in the usual way.

[statusHierarchy, thresholdValue, thresholdUncertainty] = getStatusHierarchy(network, statusObs);

nodewiseAccuracy = sum(statusHierarchy == statusObs) ./ sum(~isnan(statusObs));
timewiseAccuracy = (statusHierarchy == statusObs) * network.length ./ (~isnan(statusObs) * network.length);
totalAccuracy = sum((statusHierarchy == statusObs) * network.length) / sum(~isnan(statusObs) * network.length);


%% PLOT
plotStatus(nodeId, dates, statusObs, statusHierarchy, thresholdValue, thresholdUncertainty);
plotHierarchyGraph(hierarchyGraph)

%% FUNCTIONS
function dates = getDates(network)
	dates = network.data{1}.date;

	for n = 2:height(network)
		nodeDates = network.data{n}.date;
		newDates = nodeDates(~ismember(nodeDates, dates));

		dates = [dates; newDates]; %#ok<AGROW> 
		dates = sort(dates);
	end

end

function statusObs = getStatusMatrix(network, dates)
	% Gets data from the provided table and transforms it into matrix form.
	% One row per day and one column per node.
	statusObs = NaN(numel(dates), height(network));
	
	for s = 1:height(network)
		data = network.data{s};
	
		[id, idx] = ismember(data.date, dates);
		statusObs(idx(idx>0), s) = data.status(id);
	end
end

function O = getObservationMatrix(statusObs)
	N = size(statusObs, 2); % number of nodes
	
	O = NaN(N, N);
	for n1 = 1:N
		for n2 = 1:N
			O(n1, n2) = sum(statusObs(:, n1) == 0 & statusObs(:, n2) == 1);
		end
	end
end

function hierarchyMatrix = getHierarchyMatrix(observationMatrix)
	% net observations
	hierarchyMatrix = observationMatrix' - observationMatrix;
	hierarchyMatrix(hierarchyMatrix < 0) = 0;
end

function H = removeLoops(H)
	maxNumCycles = 10^2;
	siz = size(H);

	% remove self loops
	H = H - H';
	H(H < 0) = 0;

	% remove cycles
	for maxCycleLength = 2:2:siz(1)
		disp(strcat("Removing cycles with size ", string(maxCycleLength)));
		% dot = 0;

		while true
			Hgraph = digraph(H);
			loops = allcycles(Hgraph, 'MaxNumCycles', maxNumCycles, 'MaxCycleLength', maxCycleLength);

			% disp(loops);
			% pause();

			if isempty(loops)
				break;
			end

			H = rmLoops(H, loops, siz);

			% dot = dot + 1;
			% if dot == 100; dot = 1; fprintf('\n'); end
			% fprintf('.');
		end
	end


	function H = rmLoops(H, loops, siz)
		% for each loop
		for i = 1:numel(loops)
			% get link strengths
			loop = [loops{i} loops{i}(1)];
			
			rr = loop(1:end-1);
			cc = loop(2:end);
			ii = (cc-1)*siz(1) + rr;
			
			strength = H(ii);
			
			% remove weaker link
			% disp(' ');
			% disp(loop);
			% disp(strength);

			[~, weakerPos] = min(strength);
			H(ii(weakerPos)) = 0;
		end

		% pause();
	end
end

function [statusHie, thresholds, uncertainty] = getStatusHierarchy(network, statusObs)
	T = size(statusObs, 2); % number of nodes, i.e. possible thresholds 
	N = size(statusObs, 1); % number of surveys

	stretchLength = network.length;
	
	errorHie = NaN(N, T); % each row contains the error relative to survey N, for each possible threshold
	for threshold = 0:T
		statusTemp = [true(N, threshold), false(N, T-threshold)];
		errorHie(:, threshold+1) = (statusTemp ~= statusObs) * stretchLength;
	end

	% find treshold range associated to minimum error
	thresholdMin = NaN(N, 1);
	thresholdMax = NaN(N, 1);

	for n = 1:N
		isMinError = errorHie(n, :) == min(errorHie(n, :));

		thresholdMin(n) = find(isMinError, 1, "first") - 1;
		thresholdMax(n) = find(isMinError, 1, "last") - 1;
	end
	
	thresholds = round((thresholdMin + thresholdMax) / 2);
	uncertainty = thresholdMax - thresholdMin;

	% apply selected threshold
	statusHie = NaN(N, T);

	for n = 1:N
		statusHie(n, :) = [true(1, thresholds(n)), false(1, T-thresholds(n))];
	end
end

%% PLOT FUNCTIONS
function plotStatus(nodeId, dates, statusObs, statusHie, thresholdValue, thresholdUncertainty)
	ids = thresholdUncertainty < 20;
	dates = dates(ids);
	statusObs = statusObs(ids, :);
	statusHie = statusHie(ids, :);
	thresholdValue = thresholdValue(ids, :);

	N = size(statusHie, 2); % number of nodes
	T = size(statusHie, 1); % number of timesteps

	% preparare threshold plotting
	xHie = [thresholdValue'; thresholdValue'] + 0.5;
	yHie = [0:T-1; 1:T] + 0.5;

	% prepare confusion matrix
	C = NaN(T, N);
	C(statusObs == 1   & statusHie == 1) = 1; % true positive
	C(isnan(statusObs) & statusHie == 1) = 2; % positive, not observed
	C(statusObs == 0   & statusHie == 1) = 3; % false positive
	C(statusObs == 1   & statusHie == 0) = 4; % false negative
	C(isnan(statusObs) & statusHie == 0) = 5; % negative, not observed
	C(statusObs == 0   & statusHie == 0) = 6; % true negative

	% figure
	subplots(1, 'Margin', [10 80 20 10]);

	% plot timeseries
	imagesc(C);
	plot(xHie(:), yHie(:), '-', 'Color', color.black.rgb, 'LineWidth', 1.5);

	% x axis
	xticks(1:numel(nodeId));
	xticklabels(string(nodeId));
	xlabel('Node');
	xlim([0.5, numel(nodeId)+0.5]);

	% y axis
	% yticks(1:numel(dates));
	% yticklabels(string(dates));
	ylabel('Time');
	yticks([]);
	ylim([0.5, numel(dates)+0.5]);

	% colors and legend
	colormap([color.blue.rgb; color.lighter_blue.rgb; color.red.rgb; color.dark_blue.rgb; color.lighter_orange.rgb; color.orange.rgb]);
	cb = colorbar();
	cb.Ticks = 1+(5/6/2) :5/6: 6-(5/6/2);
	cb.TickLabels = ["TP", "P", "FP", "FN", "N", "TN"];
	cb.TickLength = 0;	
end

function plotHierarchyGraph(Hgraph)
	HreducedGraph = transreduction(Hgraph);

	figure();
	plot(HreducedGraph, 'Layout', 'layered', 'Direction', 'right', 'NodeColor', HreducedGraph.Nodes.color, 'EdgeColor', color.dark_grey.rgb, 'MarkerSize', 6);
end
