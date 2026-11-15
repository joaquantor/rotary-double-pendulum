%% Generování CSV dat z experimenty.mat
% Ujistěte se, že soubor 'experimenty.mat' je v aktuální složce
if isfile('experimenty.mat')
    load('experimenty.mat');
else
    error('Soubor experimenty.mat nebyl nalezen!');
end

% 1. Extrakce dat pro rameno (zakladna)
% V Pythonu to byl index 1, v MATLABu (indexování od 1) je to sloupec 2
% Struktura: kyvadlo2vol -> signals -> values
raw_data = kyvadlo2.signals.values(:, 2); 
time_data = kyvadlo2.time;

% 2. Přepočet na radiány
% Enkodér 500 CPR, 4x dekódování = 2000 tiků/ot
% 1000 tiků = pi radiánů
data_rad = raw_data * (pi / 1000);

% 3. Vytvoření tabulky a uložení do CSV
% Názvy sloupců musí odpovídat tomu, co čte simulační skript ('Time', 'Arm_Angle_Rad')
T = table(time_data, data_rad, 'VariableNames', {'Time', 'Arm_Angle_Rad'});

writetable(T, 'kyvadlo2_rameno_data.csv');

disp('Soubor kyvadlo_rameno_data.csv byl úspěšně vytvořen.');
disp('Nyní můžete spustit hlavní simulační skript.');