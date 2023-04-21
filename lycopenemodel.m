%original
load('allData_4Feb_2022.mat', 'model')
%model.ub(find(ismember(model.rxns,"ATPM")))=1000;
model.lb(find(ismember(model.rxns,"ATPM"))) =0;
model.ub(find(ismember(model.rxns,'EX_cbl1(e)')))=100;
[irrevmodel, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);

OR

%change original model exc rxns to reversible
excsetlb={} %upload into workspace
load('allData_4Feb_2022.mat', 'model')
model.lb(find(ismember(model.rxns,excsetlb)))=-1000;
model.ub(find(ismember(model.rxns,'EX_cbl1(e)')))=100;
%model.lb(find(ismember(model.rxns,"ATPM"))) =0;
%check exc rxns manually 
[irrevmodel, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);
% add ub=0 for _b reactions with lb=0 in original 
irrevmodel.ub(find(ismember(irrevmodel.rxns,strcat(excsetlb,'_b'))))=0;

% add reactions for lycopene pathway
irrevmodel = addReaction(irrevmodel,'FRTTs','reactionName','Farnesyltranstransferase','reactionFormula','frdp[c] + ipdp[c]  -> ppi[c] + ggdp[c]','lowerBound',0,'upperBound',1000, 'subSystem', 'Terpenoid biosynthesis (Monoterpenoid biosynthesis)')
irrevmodel = addReaction(irrevmodel,'PHYTS1s','reactionName','phytoene synthase (phdp forming)','reactionFormula','2 ggdp[c]  -> ppi[c] + phdp[c]','lowerBound',0,'upperBound',1000, 'subSystem', 'Terpenoid biosynthesis (Carotenoid biosynthesis)')
irrevmodel = addReaction(irrevmodel,'PHYTS2s','reactionName','phytoene synthase (phyt forming)','reactionFormula','phdp[c]  -> ppi[c] + phyt[c]','lowerBound',0,'upperBound',1000, 'subSystem', 'Terpenoid biosynthesis (Carotenoid biosynthesis)')
irrevmodel = addReaction(irrevmodel,'PHYTDSs','reactionName','15-cis-phytoene desaturase','reactionFormula','fad[c] + phyt[c]  -> fadh2[c] + phytfl[c]','lowerBound',0,'upperBound',1000, 'subSystem', 'Terpenoid biosynthesis (Carotenoid biosynthesis)')
irrevmodel = addReaction(irrevmodel,'PHYTFLDSs','reactionName','all-trans phytofluene desaturase','reactionFormula','fad[c] + phytfl[c]  -> fadh2[c] + z-carot[c]','lowerBound',0,'upperBound',1000, 'subSystem', 'Terpenoid biosynthesis (Carotenoid biosynthesis)')
irrevmodel = addReaction(irrevmodel,'ZCARTDSs','reactionName','all-trans-zeta-carotene desaturase ','reactionFormula','fad[c] + z-carot[c]  -> fadh2[c] + nrsprn[c]','lowerBound',0,'upperBound',1000, 'subSystem', 'Terpenoid biosynthesis (Carotenoid biosynthesis)')
irrevmodel = addReaction(irrevmodel,'CAROTDSs','reactionName','Carotene 7,8-desaturase ','reactionFormula','fad[c] + nrsprn[c]  -> fadh2[c] + lycop[c] ','lowerBound',0,'upperBound',1000, 'subSystem', 'Terpenoid biosynthesis (Carotenoid biosynthesis)')
irrevmodel = addReaction(irrevmodel,'EX_LYCOP','reactionName','Lycopene exchange','reactionFormula','lycop[c] -> ','lowerBound',0,'upperBound',1000, 'subSystem', 'Terpenoid biosynthesis (Carotenoid biosynthesis)')

% set lb for lycopene
irrevmodel.lb(2965)=0.001

% upload regulator sign file to set TR signs
data = readtable("RegulatorSign_full.csv");
regtype = data.Sign;
regulator = data.TF_KEGGID;
targets = data.GENE_KEGGID;

targets_mod = targets;
regulator_mod = regulator;
regtype_mod = regtype;
%Mapping TR target genes in metabolic model
ModGenesInTargets = targets_mod(find(ismember(targets_mod,irrevmodel.genes)));
TRsOfModGenes = regulator_mod(find(ismember(targets_mod,irrevmodel.genes)));
UniqTRsOfModGenes=unique(TRsOfModGenes);
UniqModGenesInTargets = unique(ModGenesInTargets);

RegulatorySign = regtype_mod(find(ismember(targets_mod,irrevmodel.genes)));

TargGenesInModInd = [];
for i=1:length(ModGenesInTargets)
    TargGenesInModInd = [TargGenesInModInd;find(ismember(irrevmodel.genes,ModGenesInTargets(i)))];
end

TRpseudoMetNamesEx = strcat(UniqTRsOfModGenes,'_TR[e]');
TRrxnNames = strcat(UniqTRsOfModGenes,'_transp');


% get indices of 161 TROfModGenes to group ModGenesInTargets for grouping
[~,TRgroupInd] = ismember(TRsOfModGenes,UniqTRsOfModGenes);

% group ind based on TR 
groupedModGenesInTargetsind = cell(length(UniqTRsOfModGenes),1);
for i = 1:length(UniqTRsOfModGenes)
	groupedModGenesInTargetsind{i}= TargGenesInModInd(find(TRgroupInd==i));
end
GeneName = cell(length(UniqTRsOfModGenes),1);
for k=1:length(groupedModGenesInTargetsind)
	for l=1:length(groupedModGenesInTargetsind{k})
		GeneName{k} = [GeneName{k};repelem(ModGenesInTargets(groupedModGenesInTargetsind{k}(l)),length(find(irrevmodel.rxnGeneMat(:,groupedModGenesInTargetsind{k}(l)))))'];
		end
end


% get rxn indices for genes
rxnsWithTargetGenesInds = cell(length(UniqTRsOfModGenes),1);
for k=1:length(groupedModGenesInTargetsind)
	for l=1:length(groupedModGenesInTargetsind{k})
    		rxnsWithTargetGenesInds{k} = [rxnsWithTargetGenesInds{k};find(irrevmodel.rxnGeneMat(:,groupedModGenesInTargetsind{k}(l)))];
	end
end


rxnsWithTargetGenesInds_all = cell(1,length(TargGenesInModInd));
for j=1:length(TargGenesInModInd)
    rxnsWithTargetGenesInds_all{j} = find(irrevmodel.rxnGeneMat(:,TargGenesInModInd(j)));
end

%%% get RegulatorySign for each reaction (can combine with on top)
rxnsforRegulatorySignInds = cell(length(UniqTRsOfModGenes),1);
for k=1:length(groupedModGenesInTargetsind)
	for l=1:length(groupedModGenesInTargetsind{k})
		rxnsforRegulatorySignInds{k} = [rxnsforRegulatorySignInds{k};repelem(RegulatorySign(groupedModGenesInTargetsind{k}(l)),length(find(irrevmodel.rxnGeneMat(:,groupedModGenesInTargetsind{k}(l)))))'];
		end
end

sum = cell(length(UniqTRsOfModGenes),1);
% sum regulatory sign if TR-rxn
for i = 1:length(rxnsWithTargetGenesInds)
	temp = accumarray(rxnsWithTargetGenesInds{i},rxnsforRegulatorySignInds{i});
	sum{i} = temp(rxnsWithTargetGenesInds{i})
	remove = find(sum{i}==0)
%remove if sum=0
	rxnsWithTargetGenesInds{i}(remove)=[];
	rxnsforRegulatorySignInds{i}(remove)=[];
	sum{i}(remove) = [];
end

% find rxn names from indices (for both nested loop and unnested)
rxnsWithTargetGenes = cell(length(UniqTRsOfModGenes),1);
for k=1:length(groupedModGenesInTargetsind)
	rxnsWithTargetGenes{k} = irrevmodel.rxns(rxnsWithTargetGenesInds{k},1);
end

rxnsWithTargetGenes_all = cell(length(ModGenesInTargets),1);
for k=1:length(ModGenesInTargets)
	rxnsWithTargetGenes_all{k} = irrevmodel.rxns(rxnsWithTargetGenesInds_all{k},1);
end


% add strcats 
irrevmodel.pMetNames = strcat(irrevmodel.rxns,'_R[c]');
RpseudoMetNames = cell(length(UniqTRsOfModGenes),1);
for k=1:length(groupedModGenesInTargetsind)
	RpseudoMetNames{k} = irrevmodel.pMetNames(rxnsWithTargetGenesInds{k},1);
end

uniqR = cell(length(UniqTRsOfModGenes),1);
uniqregsign = cell(length(UniqTRsOfModGenes),1);
for i = 1:length(RpseudoMetNames)
[uniqR{i},ia,ic] = unique(RpseudoMetNames{i})
uniqregsign{i} = sum{i}(ia)
end

% save point
model1 = irrevmodel;


%% add reactions (TR + R --> Z  or TR --> R)

%if <0 is repressor TR+R-> Z
%else >0 TR-> R
for m = 1:length(TRrxnNames)
    for n = 1:length(uniqR{m,1})
        TRrxnNameTemp = strcat(TRrxnNames(m),"_" ,uniqR{m,1}(n));
        TRpseudoMetNamesExTemp = TRpseudoMetNamesEx(m);
        RpseudoMetNamesTemp = uniqR{m,1}(n);
	ZpseudoMetNamesTemp = strcat("Z_" + uniqR{m,1}(n));

	if uniqregsign{m}(n)>0
		model1 = addReaction(model1, TRrxnNameTemp{1},'metaboliteList',{TRpseudoMetNamesExTemp{1},RpseudoMetNamesTemp{1}},'stoichCoeffList',[-1,1], 'reversible',false,'lowerBound',0);
	else uniqregsign{m}(n)<0
		model1 = addReaction(model1, TRrxnNameTemp{1},'metaboliteList',{TRpseudoMetNamesExTemp{1},RpseudoMetNamesTemp{1},ZpseudoMetNamesTemp{1}},'stoichCoeffList',[-1,-1,1], 'reversible',false,'lowerBound',0, 'upperBound',10);
	end
end
end



%  find ind of R in model
RmetInds = cell(161,1);
for m = 1:length(RpseudoMetNames)
    for n = 1:length(RpseudoMetNames{m,1})
		RmetInds{m}(n) = find(ismember(model1.mets,RpseudoMetNames{m,1}(n)));
    end
end

% adding R to reactions in model (A+B -> C to A+B+R -> C)
for m = 1:length(RmetInds)
    for n = 1:length(RmetInds{m,1})
		model1.S(RmetInds{m}(n),rxnsWithTargetGenesInds{m}(n)) = -1;
		end
end


% for TR exchange
TRpseudoMetNamesEx=unique(TRpseudoMetNamesEx);
% added this to name TR exchange rxns 
TRexrxnNames = strcat('EX_TR_', TRrxnNames);
% gene name for exchange need unique
TRexGeneName = cell(161,1);
for k=1:length(GeneName)
	TRexGeneName{k} = unique(GeneName{k});	
end

% adding IRREV TR exchange ( -> TR )
for l = 1:length(TRexrxnNames)
	TRexrxnNamesTemp = TRexrxnNames(l);
	TRpseudoMetNamesExTemp = TRpseudoMetNamesEx(l);
	model1 = addReaction(model1, TRexrxnNamesTemp{1},'metaboliteList',{TRpseudoMetNamesExTemp{1}},'stoichCoeffList',[1], 'reversible',false, 'lowerBound', 0, 'upperBound', 10);
end

% save point
model2=model1;

%% adding flux reactions only for repressor TRs ( -> R) and (Z->)

%find inds of -1 from rxnsforRegulatorySignInds
%arxnind = find(vertcat(rxnsforRegulatorySignInds{:})==-1);

%find inds of negative signs from rxnsforRegulatorySignInds
arxnind = find(vertcat(uniqregsign{:})<0)

%get rxnnames rxnsWithTargetGenes
%arxnname = vertcat(rxnsWithTargetGenes{:});

arxnname = vertcat(uniqR{:})
RepRpseudoMetNames=unique(arxnname(arxnind))
%need the names only for later
Repfluxrxnnames = extractBefore(RepRpseudoMetNames,'_R[c]')
% name of flux eqn for repressor TR (->R)
RdelfluxrxnNames = strcat(RepRpseudoMetNames,'_PseudoRdeltransp');
% name of sink eqn for repressor TR (R->)
%RsinkrxnNames = strcat(Repfluxrxnnames,'_PseudoRsink');
% name of R metabolite
%RepRpseudoMetNames = strcat(Repfluxrxnnames,'_R[c]');



%->R
for p = 1:length(Repfluxrxnnames)
	model2 = addReaction(model2, RdelfluxrxnNames{p,1},'metaboliteList',{RepRpseudoMetNames{p,1}},'stoichCoeffList',[1],'reversible',false,'lowerBound',0,'upperBound',1);
end

% Z->
for p = 1:length(Repfluxrxnnames)
	ZrxnNameTemp = strcat(RepRpseudoMetNames{p,1} + "_PseudoZtransport");
	ZpseudoMetNamesTemp = strcat("Z_" + RepRpseudoMetNames{p,1});
    model2 = addReaction(model2, ZrxnNameTemp{1},'metaboliteList',{ZpseudoMetNamesTemp{1}},'stoichCoeffList',[-1], 'reversible',false,'lowerBound',0,'upperBound', 10);
end

% R-> (only repressor so wrong)
%for p = 1:length(Repfluxrxnnames)
%	model2 = addReaction(model2, RsinkrxnNames{p,1},'metaboliteList',{RepRpseudoMetNames{p,1}},'stoichCoeffList',[-1],'reversible',false,'lowerBound',0,'upperBound',1000);
%end

%for all R-> (correct one is to add for all)
allRNames =strcat(unique(arxnname))
RsinkNames = strcat(allRNames,'_PseudoRsink');
for m = 1:length(allRNames)
	model2 = addReaction(model2, RsinkNames{m},'metaboliteList',{allRNames{m}},'stoichCoeffList',[-1], 'reversible',false,'lowerBound',0, 'upperBound',1000);
end

% save point


% get pfba of original model
[GeneClasses RxnClasses irrevmodelpfba] = pFBA(model);
[minfluxrev, maxfluxrev] = fluxVariability(irrevmodelpfba,'optPercentage',100,'osenseStr','max','allowLoops',true);

%%mappfba to irrev
pfba2irrev = cell(length(irrevmodelpfba.rxns)-1,1);
for i = 1:length(irrevmodelpfba.rxns)-1
	pfba2irrev{i} = find(strcmp(irrevmodel.rxns,irrevmodelpfba.rxns(i)));
end

maxfluxirrev=cell(length(irrev2rev),1); %maxfluxirrev=cell(2957,1);
minfluxirrev=cell(length(irrev2rev),1); %minfluxirrev=cell(2957,1);

tf = cellfun('isempty',maxfluxirrev) % true for empty cells
maxfluxirrev(tf) = {0}  
tf = cellfun('isempty',minfluxirrev) % true for empty cells
minfluxirrev(tf) = {0}  

pat = ("_f"|"_b")
for i =1:length(pfba2irrev)
	if contains(irrevmodelpfba.rxns(i) ,pat)
		minfluxirrev{pfba2irrev{i}}=0;
		if minfluxrev(i) <0
			maxfluxirrev{pfba2irrev{i}}= -minfluxrev(i);
		else
			maxfluxirrev{pfba2irrev{i}} = 0;
		end
		if maxfluxrev(i) <0
			maxfluxirrev{pfba2irrev{i}} = 0;
		else 
			maxfluxirrev{pfba2irrev{i}} = maxfluxrev(i);
		end

	else 
		minfluxirrev{pfba2irrev{i}}=0;
		maxfluxirrev{pfba2irrev{i}}=maxfluxrev(i);	

	end
end

% repressor rxns ind from modelrxns
corr = cell(length(Repfluxrxnnames),1);
for i = 1:length(Repfluxrxnnames)
	corr{i} = find(strcmp(irrevmodelpfba.rxns,Repfluxrxnnames(i)));
end

%% set fluxes for irrev based on rev

[originalminfluxrev, originalmaxfluxrev] = fluxVariability(model,'optPercentage',100,'osenseStr','max','allowLoops',true);
originalmaxfluxirrev=cell(2957,1);
originalminfluxirrev=cell(2957,1);
for i =1:length(rev2irrev)
	if length(rev2irrev{i})==1
		originalminfluxirrev{rev2irrev{i}}=0;
		originalmaxfluxirrev{rev2irrev{i}}=originalmaxfluxrev(i);
	else 
		tempfind = rev2irrev{i}(1);
		tempbind = rev2irrev{i}(2);
		originalminfluxirrev{tempfind}=0;
		originalminfluxirrev{tempbind}=0;
		if originalmaxfluxrev(i) <0
			originalmaxfluxirrev{tempfind} = 0;
		else 
			originalmaxfluxirrev{tempfind} = originalmaxfluxrev(i);
		end
		if originalminfluxrev(i) <0
			originalmaxfluxirrev{tempbind}= -originalminfluxrev(i);
		else
			originalmaxfluxirrev{tempbind} = 0;
		end
	end
end

% process results and get high flux
originalcorr = cell(length(Repfluxrxnnames),1);
for i = 1:length(Repfluxrxnnames)
	originalcorr{i} = find(strcmp(irrevmodel.rxns,Repfluxrxnnames(i)));
end

listoriginal = zeros(1,length(originalcorr));
for i=1:length(originalcorr)
listoriginal(i) = originalmaxfluxirrev{originalcorr{i}};
end

highfluxind = find(listoriginal > 100);
%highfluxnames = Repfluxrxnnames(find(listoriginal > 100))


% change ub to 10% flux if >200 for ->R  1114 only 
% change lb to 90% flux if >200 for Z->  1114 only 

model3=model2;
for i= 1:length(originalcorr)
if isempty(originalcorr{i})
continue
else
	if (originalmaxfluxirrev{originalcorr{i}})==0
		model3.ub(length(irrevmodel.rxns)+length(vertcat(uniqR{:}))+length(TRexrxnNames)+i) = 0;
		model3.lb(length(irrevmodel.rxns)+length(vertcat(uniqR{:}))+length(TRexrxnNames)+length(Repfluxrxnnames)+i) = 0;
	elseif (originalmaxfluxirrev{originalcorr{i}})<10
		model3.ub(length(irrevmodel.rxns)+length(vertcat(uniqR{:}))+length(TRexrxnNames)+i) = originalmaxfluxirrev{originalcorr{i}};
		model3.lb(length(irrevmodel.rxns)+length(vertcat(uniqR{:}))+length(TRexrxnNames)+length(Repfluxrxnnames)+i) = 0.1*originalmaxfluxirrev{originalcorr{i}};
	elseif (originalmaxfluxirrev{originalcorr{i}})<100
		model3.ub(length(irrevmodel.rxns)+length(vertcat(uniqR{:}))+length(TRexrxnNames)+i) = 0.1*originalmaxfluxirrev{originalcorr{i}};
		model3.lb(length(irrevmodel.rxns)+length(vertcat(uniqR{:}))+length(TRexrxnNames)+length(Repfluxrxnnames)+i) = 0.01*originalmaxfluxirrev{originalcorr{i}};
	else
		model3.ub(length(irrevmodel.rxns)+length(vertcat(uniqR{:}))+length(TRexrxnNames)+i) = originalmaxfluxirrev{originalcorr{i}};
		model3.lb(length(irrevmodel.rxns)+length(vertcat(uniqR{:}))+length(TRexrxnNames)+length(Repfluxrxnnames)+i) = 0;
	end
end
end
% if highflux need to separately fix
for ii=1:length(highfluxind)
	model3.ub(length(irrevmodel.rxns)+length(vertcat(uniqR{:}))+length(TRexrxnNames)+highfluxind(ii)) = maxfluxirrev{corr{highfluxind(ii)}}
	model3.lb(length(irrevmodel.rxns)+length(vertcat(uniqR{:}))+length(TRexrxnNames)+length(Repfluxrxnnames)+highfluxind(ii)) = 0.1*maxfluxirrev{corr{highfluxind(ii)}};
end


% save point for model3



%% relaxation (manually check in excel)
load('relax_rxns.mat')

model3.lb(find(ismember(model3.rxns,"EX_LYCOP")))=0; %model3.lb(2965)=0
model3.lb(find(ismember(model3.rxns,"ATPM")))=8.39; %model3.lb(375)= 8.39
%find maxed rxns (

modelt=model3
modelt.ub(find(ismember(modelt.rxns,'EX_TR_b0683_transp')))=100; 
modelt.ub(find(ismember(modelt.rxns,'EX_TR_b0889_transp')))=100; 
modelt.ub(find(ismember(modelt.rxns,'EX_TR_b1221_transp')))=100;
modelt.ub(find(ismember(modelt.rxns,rxnst)))=10*ones(57,1); %rxnst is manually found from filtering in excel
optimizeCbModel(modelt)

model4=modelt;
model4.ub(find(ismember(model4.rxns,'EX_glc(e)_b')))=8.5;
model4.ub(find(ismember(model4.rxns,'EX_o2(e)_b')))=14.6;
optimizeCbModel(model4)

% second relaxation
% don't relax b0080

modelt2=model4
modelt2.ub(find(ismember(modelt2.rxns,rxnst2)))=10*ones(18,1);  %rxnst2 is manually found from filtering in excel
modelt2.ub(find(ismember(modelt2.rxns,'EX_TR_b1827_transp')))=100; 
modelt2.ub(find(ismember(modelt2.rxns,'EX_TR_b3357_transp')))=100;
optimizeCbModel(modelt2)

% third relaxation

modelt4=modelt2
%modelt4.ub(find(ismember(modelt4.rxns,'EX_TR_b4062_transp')))=100;
modelt4.ub(find(ismember(modelt4.rxns,rxnst3)))=100*ones(7,1);  %rxnst3 is manually found from filtering in excel
modelt4.ub(find(ismember(modelt4.rxns,'UGLYCH_R[c]_PseudoRdeltransp')))=10*ones(1,1);
optimizeCbModel(modelt4)

%CHANGE TO RXN NAME 
modelt4.ub(find(ismember(modelt4.rxns,"b0683_transp_GLUDy_b_R[c]")))=1;
modelt4.ub(find(ismember(modelt4.rxns,"b1988_transp_GLUDy_b_R[c]")))=1;
modelt4.ub(find(ismember(modelt4.rxns,"b2916_transp_GLUDy_b_R[c]")))=1;
modelt4.ub(find(ismember(modelt4.rxns,"b3357_transp_GLUDy_b_R[c]")))=1;
optimizeCbModel(modelt4)

modelt4.ub(find(ismember(modelt4.rxns,"CYTBD2pp_R[c]_PseudoRsink")))=0
modelt4.ub(find(ismember(modelt4.rxns,"CYTBDpp_R[c]_PseudoRsink")))=1000
modelt4.ub(find(ismember(modelt4.rxns,"HYD1pp_R[c]_PseudoRsink")))=0
modelt4.ub(find(ismember(modelt4.rxns,"HYD2pp_R[c]_PseudoRsink")))=0
modelt4.ub(find(ismember(modelt4.rxns,"HYD3pp_R[c]_PseudoRsink")))=0
modelt4.ub(find(ismember(modelt4.rxns,"NTP3pp_R[c]_PseudoRsink")))=0
modelt4.ub(find(ismember(modelt4.rxns,"PHYTSpp_R[c]_PseudoRsink")))=0

% change from using CYTBDpp rxn to using CYTBO3_4pp rxn
modelt4.S(find(ismember(modelt4.mets,"CYTBDpp_R[c]")),find(ismember(modelt4.rxns,"CYTBDpp")))=0; %modelt4.S(1757,516)=0;
modelt4.S(find(ismember(modelt4.mets,"CYTBDpp_R[c]")),find(ismember(modelt4.rxns,"CYTBO3_4pp")))=-1; %modelt4.S(1757,517)=-1;

% set to 100 bc the flux for the other rxn is ~29
modelt4.ub(find(ismember(modelt4.rxns,'EX_TR_b0564_transp')))=100; 
optimizeCbModel(modelt4)


% save point for modelt4


% modelt4 0.6637 (427 rxns)
excrxn1={} % upload to workspace

% to see lycopene flux 
perm100flux = cell(1000,2);
modelx=modelt4;
modelx.lb(find(modelx.c))=0.6337*0.5
modelx.c(find(modelx.c))=0;
modelx.c(find(ismember(modelx.rxns,'EX_LYCOP')))=1;
modelx.ub(find(ismember(modelx.rxns,'EX_TR_b0080_transp')))=100;
for ii= 1:length(excrxn1)
	modelx.ub(find(ismember(modelx.rxns,excrxn1(ii))))=0;
	optimizeCbModel(modelx)
	perm100flux(i,1)= num2cell(ans.f);
	perm100flux(i,2)=num2cell(ans.x(1005));
end
