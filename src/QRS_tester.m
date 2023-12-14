%% Funkce pro testování úspěšnosti detekce QRS
% Vstupy: 
% Ref - referenční pozice z databáze
% Poz - Vámi detekované pozice
% fs - vzorkovací frekvence signálu
% Výstupy:
% ACC - celková úspěšnost detekce
% Se - senzitivita
% PP - Pozitivní prediktivita
% TP - true positive
% FP - false positive
% FN - false negative
%%

function [ACC, Se, PP, TP, FP, FN] = QRS_tester(ref, poz, fs)

tol_ms =100;  
tol_vz = round(tol_ms/1000*fs);  

    [TP,FP,FN,dif, Fpp, Fnn] = SePms(poz,ref,tol_vz);
    
    
    Se = (TP/(TP+FN))*100;
    PP = (TP/(TP+FP))*100;
    ACC=(Se+PP)/2;

    m = mean(dif)/fs*1000;   % ms
    s = std(dif)/fs*1000;    % ms
Fpp=Fpp;
Fnn=Fnn;
end

%------ POMOCNÉ FUNKCE ----------------------------------------------------

function [TP,FP,FN,dif, Fpp, Fnn]= SePms(QRS,refQRS,tol)
% clear all; close; clc
% refQRS =    [100 200 300 400 420];
% QRS =       [85 215 250 301];
% tol = 30;

TP = 0;
FP = 0;
Fpp=[];
Fnn=[];
FN = 0;
dif = 0;  
j = 0;

for i = 1:length(QRS)
    pom = find(QRS(i)>=refQRS-tol & QRS(i)<=refQRS+tol, 1);
    if isempty(pom)
        FP = FP+1;
        Fpp=[Fpp i];
    end
end

for i = 1:length(refQRS)
    pom = find(refQRS(i)>=QRS-tol & refQRS(i)<=QRS+tol);
    if isempty(pom)
        FN = FN+1;
                Fnn=[Fnn i];

    elseif  ~isempty(pom)
        TP = TP+1;
        j = j+1;
        idx = find( abs(refQRS(i)-QRS(pom)) == min(abs(refQRS(i)-QRS(pom))), 1, 'first');
        dif(j) = QRS(pom(idx)) - refQRS(i);
        if length(pom) >= 2
            FP = FP + length(pom) - 1;
                    Fpp=[Fpp i];

        end
    end
end

if isempty(refQRS)
    TP = 0;
    FP = 0;
    FN = 0;
    dif = 0;
end
end