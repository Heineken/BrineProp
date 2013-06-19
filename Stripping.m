path(path,'d:\Program Files (x86)\Dymola 2013\Mfiles')
path(path,'d:\Program Files (x86)\Dymola 2013\Mfiles\dymtools')

%%pick file and time
[filename,pathname] = uigetfile('..\*.mat','Select the M-file'); %choose mat-file
clear('-regexp', '?!(^filename|^pathname)'); %if not cancelled clear all except ...
if ~filename, return, end
res = dymload(fullfile(pathname, filename)); %load mat-file

%% extract
n=[1 0]*dymget(res, 'n');
X_g_cell=dymget(res, 'X_g');
X_g=zeros(n,4);
p_cell=dymget(res, 'p');
p=zeros(n,1);
GVF=zeros(n,1);
GVF=zeros(n,9);
for i=1:n
    p(i)=p_cell{i}(1);
    GVF(i)=[1 0]*dymget(res, ['props[' num2str(i) '].state.GVF']);
    X(i,:)=dymget(res, ['props[' num2str(i) '].X']);
    for j=6:size(X_g_cell,2)
        X_g(i,j-5)=X_g_cell{i,j}(1);
    end
end

%% Plot
plot([X_g GVF X_g(:,2)./X_g(:,3)], p)
set(gca,'YDir', 'reverse');

legend('CO2', 'N2', 'CH4','H2O','GVF')