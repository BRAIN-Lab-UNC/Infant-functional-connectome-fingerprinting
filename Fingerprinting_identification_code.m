% Code for feature-selection based identification test
% Citation: Dan Hu, Fan Wang, Han Zhang, Zhengwang Wu, Zhen Zhou, Guoshi Li, Li Wang, Weili Lin, Gang Li. Existence of Functional Connectome Fingerprint During Infancy and Its Stability Over Months. Journal of Neuroscience 2022, 42(3): 377-389
% Contact: danhu0055@gmail.com, gang_li@med.unc.edu
Data1=Base_Set;
Data2=Target_Set;
RepTime=20;
K=10;
t=1;
for i=1:RepTime
    indices = crossvalind('Kfold',size(Data1,1),K);
    for ii=1:K
        TrainData1=Data1(indices~=ii,:);
        TrainData2=Data2(indices~=ii,:);
        TestData1=Data1(indices==ii,:);
        TestData2=Data2(indices==ii,:);
        TheMeasure=(std(TrainData1)+std(TrainData2))/2; 
        [FinalAlpha1,FinalAlpha2]=InnerAlpha_choose(TheMeasure,TrainData1,TrainData2);
        UsefulID=find(TheMeasure>quantile(TheMeasure,FinalAlpha1));
        TheChosen(UsefulID,t)=TheMeasure(UsefulID);
        t=t+1;
        [TestAcc1(i,ii),testLen1(i,ii)]=CalIDrate(TheMeasure,FinalAlpha1,Data1,TestData2,find(indices==ii));      
        [TestAcc2(i,ii),testLen2(i,ii)]=CalIDrate(TheMeasure,FinalAlpha2,Data2,TestData1,find(indices==ii));      
    end
    FinalTestAcc_r1(i)=sum(testLen1(i,:))/size(Data1,1)
    FinalTestAcc_r2(i)=sum(testLen2(i,:))/size(Data2,1)    
end


function [FinalAlpha1,FinalAlpha2]=InnerAlpha_choose(TheMeasure,TrainData1,TrainData2)

TestAcc1_old=0;
TestAcc2_old=0;
TheAlpha=[0.8:0.02:0.98,0.982:0.002:0.999];

for kk=1:length(TheAlpha)
    alpha=TheAlpha(kk);
    [TestAcc1,~]=CalIDrate(TheMeasure,alpha,TrainData1,TrainData2,[1:1:size(TrainData2,1)]'); 
    if TestAcc1>TestAcc1_old
        FinalAlpha1=alpha;
        TestAcc1_old=TestAcc1;
    end
    
    [TestAcc2,~]=CalIDrate(TheMeasure,alpha,TrainData2,TrainData1,[1:1:size(TrainData2,1)]');
    if TestAcc2>TestAcc2_old
        FinalAlpha2=alpha;
        TestAcc2_old=TestAcc2;
    end
end
end

function [Acc1,correct1]=CalIDrate(TheMeasure,alpha,Data1,TrainData2,testID)
    UsefulID=find(TheMeasure>quantile(TheMeasure,alpha));
    Origin_data=Data1(:,UsefulID);
    Target_data=TrainData2(:,UsefulID);
    Corr=[];
    for i=1:size(Target_data,1)

        for j=1:size(Origin_data,1)
            Temp=corrcoef((Origin_data(j,1:end))',(Target_data(i,1:end))');
            Corr(i,j)=Temp(1,2);

        end
    end
    [~,BB]=max(Corr');
    Tar=testID';
    Acc1=length(find(BB-Tar==0))/size(Target_data,1);
    correct1=length(find(BB-Tar==0));
end

