
% clear
% Data_input=input


av=mean(Data_input);
stnd=std(Data_input);
L=size(Data_input,1);
tst = icdf('normal',1-1/(4*L),0,1);
tst_data=(Data_input-av)/stnd;
[RR,C]=find(tst_data>tst);
excluded=Data_input(RR,:);
Data_input(RR,:)=[];

c19 = 1;

while isempty(excluded)~=1
    av=mean(Data_input); stnd=std(Data_input);
    L=size(Data_input,1);
    tst = icdf('normal',1-1/(4*L),0,1);
    tst_data=(Data_input-av)/stnd;
    [RR,C]=find(tst_data>tst);
    excluded=Data_input(RR,:);
    Data_input(RR,:)=[];
    c19 = 1+c19;
end
% RR