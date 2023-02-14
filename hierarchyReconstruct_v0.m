clear;
close all;
clc;
% load('viterbo_all.mat');
load('valfredda.mat');

%% restructure data and create graph

% definisco le date di cui usare le osservazioni. C'è la funzione che lo fa in base ai dati, oppure si può fare a mano.
dates = getDates(network);
% dates = datetime(2020, 01, 01):datetime(2020, 01, 12);

nodeId = network.Name;

% matrice di stati osservati: una riga per ogni giorno, una colonna per ogni nodo
statusObs = getStatusMatrix(network, dates);

% trasformo le osservazioni singole in un grafo, in cui il nodo A è connesso al nodo B se c'è un'osservazione in cui
% abbiamo visto A acceso e B spento (che ci direbbe che A è più persistente di B). Il peso dell'arco che unisce A a B è
% il numero di osservazioni di A acceso e B spento. Se una (o più) volte è stato visto A spento e B acceso, ci sarà
% anche un arco diretto da B ad A. Questo grafo è rappresentato da una matrice con N (numero total di nodi) righe ed N
% colonne, dove nella cella (r,c), dove r è la riga e c la colonna, c'è il numero di osservazioni in cui il nodo r era
% spento e c acceso.
observationMatrix = getObservationMatrix(statusObs);

% creo un secondo grafo, che descriverà la gerarchia: se abbiamo visto 5 volte A acceso e B spento, e 2 volte A spento e
% B acceso, in questo grafo ci sarà un collegamento da A a B con peso 3 (cioè il numero netto di volte in cui le
% osservazoini dicono che A è più persistente di B). Lo stesso tipo di problema può accadere anche con 3 o più nodi: le
% osservazioni ci dicono che A è più persistente di B (A > B), B > C, ma anche C > A. Questo può capitare quando le
% diverse coppie di nodi sono state osservate in tempi diversi, e crea dei loop chiusi nel grafo, che danno fastidio.
% Per questo motivo, identifichiamo tutti i loop e li apriamo eliminando l'arco che pesa meno. In questo modo resta il
% grafo che identifica la gerarchia andando contro meno osservazioni possibili.
hierarchyMatrix = getHierarchyMatrix(observationMatrix);
hierarchyMatrix = removeLoops(hierarchyMatrix);
hierarchyGraph = digraph(hierarchyMatrix, network);

% [TODO] in caso ci siano troppi loop, usare il parametro MinCycleLength: con un ciclo, partire da loop di lunghezza grande e
% toglierli, e man mano togliere i loop via via più piccoli fino ad 1.

%% reorder graph topologically
% riordino i nodi in modo topologico: in questo modo i nodi più persistenti arrivano prima di quelli meno persistenti.
% Per questa operazione serve che non ci siano loop nel grafo, per questo li abbiamo tolti.
% Nota: ci possono essere molteplici ordinamenti di nodi che soddisfano tutti i link di gerarchia; a noi ne basta uno
% qualsiasi perché i nodi che possono scambiarsi posizione tra loro sono quelli che non hanno osservazioni in comune, e
% quindi non incideranno sul calcolo dell'accuracy.
% Penso che da questa wiki sia tutto chiaro: https://en.wikipedia.org/wiki/Topological_sorting
order = toposort(hierarchyGraph);

network = network(order, :);
nodeId = nodeId(order);
statusObs = statusObs(:, order);
observationMatrix = observationMatrix(order, order);
hierarchyMatrix = hierarchyMatrix(order, order);

% ora tutti i nodi sono ordinati secondo la gerarchia.

%% check hierarchy accuracy
% possibilità 1: calcolare la frazione di osservazioni (intese come osservazione di una coppia di nodi, di cui uno
% spento ed uno acceso) che sono in accordo con la gerarchia. Questo coincide con sommare gli elementi del triangolo
% superiore di observation matrix diviso la somma totale degli elementi (una volta riordinata la matrice, tutte le
% osservazoini che seguono la gerarchia stanno nel triangolo superiore, perché i nodi più persistenti stanno prima
% nell'ordine gerarchico).

% Questo per il paper non lo guardiamo
hierarchyAccuracy = sum(tril(observationMatrix), 'all') / sum(observationMatrix, 'all');

% possibilità 2: per ogni data (ovvero ogni riga di statusObs), prendere tutte le osservazioni disponibli in statusObs
% (intese come osservazione di un singolo nodo che può essere acceso o spento). Secondo la gerarchia, i nodi si
% accendono dal primo in avanti fino ad un certo punto; tutti i nodi dopo saranno spenti. Sarà da far scorrere il numero
% di nodi accesi da 1 ad N, e per ognuno calcolare la lunghezza "true positive o negative" (per cui lo stato previsto
% dall'accensione gerarchica è lo stesso dello stato osservato), e trovare il numero di nodi accesi che massimizza
% questa true length (per la data corrente). Questo simula il trovare la soglia di persistenza che delimita i nodi
% accesi e spenti, solo che non conosciamo effettivamente le persistenze ma solo la loro gerarchia. In questo modo
% possiamo costruire una matrice statusHierarchy, analoga a statusObs ma che contiene lo stato previsto dall'accensione
% gerarchica.
% NOTA: questo metodo vale solo per calcolare l'accuracy, non per ricostruire lo stato di nodi non osservati in base a
% quelli osservati. Questo perché l'ordinamento topologico non è univoco nei casi in cui la gerarchia non sia stata
% osservata.

%    [A]      Questi nodi possono essere ordinati come ABCDFE oppure AFBCDE (o anche altri modi). Guardando solo
%   /  \	  il primo ordinamento, nel caso in cui F sia osservato acceso verrebbe da dare come accesi anche BCD, 
% [B]  [F]    cosa che però la gerarchia non garantisce (e che non sarebbe vera usando il secondo ordinamento).
%  |    |     Questo non è un problema nel calcolo dell'accuracy, perché quando abbiamo osservato F non abbiamo
% [C]   |     osservazioni di BCD (altrimenti esisterebbe una gerarchia tra F e BDC), per cui quando F conta, non
%  |    |     contano BCD, e viceversa.
% [D]   |
%   \  /
%   [E]
[statusHierarchy, thresholdValue, thresholdUncertainty] = getStatusHierarchy(network, statusObs);

% Da qui possiamo calcolare 3 accuracy diverse: 
% - node-wise accuracy: calcolata per ogni % nodo come numero di volte in cui lo stato gerarchico coincide con quello 
%   osservato, diviso il numero di osservazioni  per quel nodo 
% - time-wise accuracy: calcolata per ogni data come numero di nodi in cui lo stato gerarchico coincide con quello 
%   osservato, diviso il numero di nodi osservati in quella data 
% - total accuracy: calcolata come numero totale di osservazioni giuste  (mettendo insieme tutti i nodi e tutte le 
%   date) diviso il numero totale di osservazioni.

% questi sono validii se ogni nodo ha la stessa lunghezza
% nodewiseAccuracy = sum(statusHierarchy == statusObs) ./ sum(~isnan(statusObs));
% timewiseAccuracy = sum(statusHierarchy == statusObs, 2) ./ sum(~isnan(statusObs), 2);
% totalAccuracy = sum(statusHierarchy == statusObs, 'all') ./ sum(~isnan(statusObs), 'all');

% in generale
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
	% prende i dati dalla mega tabella e li trasforma in matrice
	% la matrice avrà una riga per ogni data, ed una colonna per ogni nodo (stazione)
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
	T = size(statusObs, 2); % numero di nodi, i.e. di possibili soglie
	N = size(statusObs, 1); % numero di survey

	stretchLength = network.length;
	
	errorHie = NaN(N, T); % ogni riga conterrà l'errore per il survey N, al variare della soglia
	for threshold = 0:T
		statusTemp = [true(N, threshold), false(N, T-threshold)];
		errorHie(:, threshold+1) = (statusTemp ~= statusObs) * stretchLength;
	end

	% trova il range di threshold associato al minimo errore
	thresholdMin = NaN(N, 1);
	thresholdMax = NaN(N, 1);

	for n = 1:N
		isMinError = errorHie(n, :) == min(errorHie(n, :));

		thresholdMin(n) = find(isMinError, 1, "first") - 1;
		thresholdMax(n) = find(isMinError, 1, "last") - 1;
	end
	
	thresholds = round((thresholdMin + thresholdMax) / 2);
	uncertainty = thresholdMax - thresholdMin;

	% ricalcola lo stato dei nodi con il threshold finale
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