%% CALCULATE FLUX SLOPE AND GET TOP EXCHANGE RXNS
% permutate 1000 times and get avg of positive slopes

% upload matrix of bounds
perm100flux = cell(1000,2); % create empty flux
q=cell(length(modelx.rxns),1000); % modelx obtained from lycopenemodel.mat code
perm100=readmatrix("randomperm100_427.csv"); % matrix of randomly generated upper bounds for the 427 reactions (0/10/100)
perm100(1,:)=[];
perm100(:,1)=[];

% get fluxslopes permuted 1000 times
list1=[];
for i = 1:1000   
 modelx=modelt4; %modelt4 obtained from lycopenemodel.mat code
 % in this case don't need but use if want to set objective, lb,ub 
 %modelx.lb(find(ismember(modelx.rxns,'Ec_biomass_iAF1260_core_59p81M')))=0.6337*0.5;
 %modelx.c(find(ismember(modelx.rxns,'Ec_biomass_iAF1260_core_59p81M')))=0;
 %modelx.c(find(ismember(modelx.rxns,'EX_LYCOP')))=1;
 %modelx.ub(find(ismember(modelx.rxns,'EX_TR_b0080_transp')))=100;
for ii= 1:length(excrxn) %excrxn can be uploaded to workspace
 if perm100(i,ii) ==0
  modelx.ub(find(ismember(modelx.rxns,excrxn(ii))))=0;
  elseif perm100(i,ii) ==10
  modelx.ub(find(ismember(modelx.rxns,excrxn(ii))))=10;
 else perm100(i,ii) ==100
  modelx.ub(find(ismember(modelx.rxns,excrxn(ii))))=100;
 end
end
sol = optimizeCbModel(modelx)
perm100flux(i,1)= num2cell(sol.f);
perm100flux(i,2)=num2cell(sol.x(1005));
 try
    q(:,i)=num2cell(calculateFluxSlope(modelx,'EX_LYCOP'));
  catch ME
    disp(i);
    list1=[list1;i]
    continue;  % Jump to next iteration of: for i
 end
end

%% process results
indexc=find(ismember(modelx.rxns,excrxn));
qexc = q(indexc,:);
qmat=cell2mat(qexc);

% select only positive values and get average of positive values
posValue=cell(length(excrxn),1);
avgValue=cell(length(excrxn),1);
for i =1:length(excrxn)
 temp=qmat(i,:);
 posValue(i,1)={temp(qmat(i,:)>0)};
 avgValue(i,1) = {mean(posValue{i})};
end

% get threshold
avgValueMatrix = cell2mat(avgValue)
avgValueMatrix(isnan(avgValueMatrix))=[];
threshold=mean(avgValueMatrix) % get threshold 

% find bp above threshold
bpValue=cell(length(excrxn),1);
for i =1:length(excrxn)
 temp=qmat(i,:)';
 bpValue(i,1)={temp(qmat(i,:)>threshold)};
end

bpind= find((~cellfun('isempty',bpValue)));
bpValuenew = bpValue(~cellfun('isempty',bpValue));
boxplotGroup(bpValuenew') % plot boxplot

% find bp above 10
bpValue_s=cell(length(excrxn),1);
for i =1:length(excrxn)
 temp=qmat(i,:)';
 bpValue_s(i,1)={temp(qmat(i,:)>10)};
end

%bpind_s = find((~cellfun('isempty',bpValue_s)));
%bpValuenew_s = bpValue_s(~cellfun('isempty',bpValue_s));
%boxplotGroup(bpValuenew_s') % plot boxplot

% see rxns from excrxn(bpind_s)
excrxn(bpind_s)


% find most number of bpValue
bpList=[];
for x = 1:length(posValue)
   bpList = [bpList,length(posValue{x})]
end
bpList_100 = find(bpList>100) %find those with >100 positive values
excrxn(bpList_100)
