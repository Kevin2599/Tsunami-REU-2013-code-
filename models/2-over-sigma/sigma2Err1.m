
maxNLambda = 80;
% measure: max run up/down
% time of run up/down
% integral of difference in eta
options = {'model_type','2/sigma', ...
	'maxl',20,'dlambda',.01,'dsigma',1,'maxsigma',15,'timeFixEnd',1, ...
	'plotTime',false,'plotLambda',false,'plotModelProgress',false,'logProgress',false};

results = cell();
for nLambda = 1:maxNLambda-1
	malprintf('\r nLambda: %d/%d',nLambda,maxNLambda);
	results{nLambda} = phiModel(options{:},'num_lambda',nLambda);
end

println('');