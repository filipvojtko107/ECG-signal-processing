clear all;

% Vektor k souborum s daty
soubory = ["101.mat", "103.mat", "106.mat", "117.mat", "119.mat", "122.mat", "214.mat", "223.mat", "231.mat", "232.mat"];
VZORKU = 650000;
VZORKOVACI_FREKVENCE = 360;
MIN_FREKV = 18;  % Minimalni hodnota frekvencni slozky signalu (high-pass)
MAX_FREKV = 35;  % Maximalni hodnota frekvencni slozky signalu (low-pass)
vysledky = zeros(1, length(soubory));
tf = zeros(1, length(soubory));

for i=1:1:length(soubory)
    % Nacteni dat signalu
    data = load(soubory(i));
    ekg_signal = data.x(1:end, 1);

    % Filtrace signalu
    ekg_signal_f = filtr(ekg_signal, VZORKOVACI_FREKVENCE, MIN_FREKV, MAX_FREKV);

    %{
    hold on;
    plot(ekg_signal);
    plot(ekg_signal_f);
    legend('Behem filtrace', 'Pred filtraci', 'Po filtraci');
    %}

    % Rozkmitani signalu
    ekg_signal_ok = add_high_frequency(ekg_signal_f);

    % Nalezeni pozic QRS komplexu
    qrs_pozice = vypocitej(ekg_signal_ok, VZORKU);

    % Vykresleni test
    %plot(ekg_signal);
    %hold on;
    %stem(qrs_pozice, ones(length(qrs_pozice), 1));

    % Provedeni testu a ulozeni vysldku testu
    [ACC, Se, PP, TP, FP, FN] = QRS_tester(data.ann, qrs_pozice, VZORKOVACI_FREKVENCE);
    vysledky(i) = ACC;
    tf(i) = length(qrs_pozice) / (VZORKU / 360 / 60);
end
vypis_vysledky(vysledky, soubory, tf);

 
% Vypis vysledku
function [] = vypis_vysledky(vysledky, soubory, tf)
    fprintf("Vysledky:\n");
    for i=1:1:length(soubory)
        fprintf("%s -> %.2f %%; Tepova frekv: %.2f\n", soubory(i), vysledky(i), tf(i));
    end
    vysledky_prumer = mean(vysledky);
    fprintf("Prumer: %.2f %%\n", vysledky_prumer);
end


% Filtrovani signalu, kombinace low-pass a high-pass filtru (band-pass filtr)
function [ekg_signal_filtr] = filtr(ekg_signal, vzorkovaci_frekvence, high_pass_meze, low_pass_meze)
    % Linearni filtrovani
    low = high_pass_meze / (vzorkovaci_frekvence / 2);
    high = low_pass_meze / (vzorkovaci_frekvence / 2);
    [b, a] = fir1(27, [low, high], "bandpass");
    ekg_signal_filtr = filtfilt(b, a, ekg_signal);

    % Zvetseni magnitudy signalu pro lepsi rozkmitani signalu
    ekg_signal_filtr = ekg_signal_filtr + 0.000001;

    % Zdurazneni signalu pomoci nelinearniho filtrovani
    for i=1:1:length(ekg_signal_filtr)
        value = sign(ekg_signal_filtr(i)) * (ekg_signal_filtr(i)^2);
        ekg_signal_filtr(i) = value;
    end
end


% Pridani vysokofrekvencniho a nizkoamplitudoveho signalu
function [ekg_signal_new] = add_high_frequency(ekg_signal)
    n = length(ekg_signal);
    forgetting_factor = 0.5;
    c = 4;
    K = zeros(1, n);
    K(1) = ekg_signal(1);
    ekg_signal_new = zeros(1, n);
    avg = mean(abs(ekg_signal));

    for i=length(ekg_signal):-1:2
        k = (forgetting_factor*K(i - 1)) + ((1 - forgetting_factor)*abs(ekg_signal(i))*c);
        K(i) = k;
    end

    for i=1:1:length(ekg_signal_new)
        K(i) = ((-1)^i) * K(i);
        if abs(ekg_signal(i)) < avg
            ekg_signal_new(i) = K(i);
        else
            ekg_signal_new(i) = ekg_signal(i);
        end
    end
end


function [qrs_pozice] = vypocitej(ekg_signal, vzorku)
    segment_size = 160;
    count = ceil(vzorku / segment_size);
    D = zeros(1, count);
    D(1) = segment_size;

    di = 2;  % D index
    si = 1;  % segment index
    for i=1:1:count
        ss = si + segment_size - 1;
        segment = [];
        if ss > vzorku
            segment = ekg_signal(ss-segment_size:end);
        else
            segment = ekg_signal(si:ss);
        end

        d = 0;  % pocet pruchodu nulou v jednom segmentu
        for n=2:1:length(segment)
            value = abs((sign(segment(n)) - sign(segment(n-1))) / 2);
            d = d + value;
        end

        D(di) = d;
        di = di + 1;
        si = si + segment_size;
    end

    % Vyjadri QRS pozice
    forgetting_factor = 0.7;
    qrs_pozice = [];
    threshold = zeros(1, count);
    threshold(1) = D(1);
    si = 1;

    for i=2:1:count
        ss = si + segment_size - 1;
        segment = [];
        if ss > vzorku
            segment = ekg_signal(ss-segment_size:end);
        else
            segment = ekg_signal(si:ss);
        end
        
        threshold(i) = forgetting_factor*threshold(i-1)+(1-forgetting_factor)*D(i);
        if D(i) < threshold(i)
           [mv, mi] = max(abs(segment));
           mi = mi + si - 1;
           qrs_pozice = [qrs_pozice, mi];
        end
        si = si + segment_size;
    end
end


