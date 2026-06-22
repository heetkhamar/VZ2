%% 1) Import-Optionen definieren
%    Automatisches erkennen der Struktur
opts = detectImportOptions('KHP15_0,6.txt','Delimiter',';'); % WICHTIG: Datensatz muss gleich wie unten sein!!!!

%    Jetzt stellen wir sicher, dass die Spalte "Legierung" als Text (char) gelesen wird.
%    Der korrekte Variablenname in 'opts' kann je nach Header etwas anders heißen.
%    Wenn in opts.VariableNames die Spalte "Legierung" heißt, dann:
opts = setvaropts(opts,'Legierung','Type','char');

%% 2) Tabelle einlesen
T = readtable('KHP15_0,6.txt',opts); % WICHTIG: Datensatz muss gleich wie oben sein!!!!

% Kontrolle:
head(T)

% 3) Spalte 'Legierung' in ein Cell-Array umwandeln
T.Legierung = cellstr(T.Legierung);

% 4) Mapping-Dictionary erstellen
legMap = containers.Map( ...
    {'CU-ETP','CU-PHC','CUSN4','CUSN5','CUSN6','CUSN6HP','CUSN8','CUSN8HP','CUSN10','CUSN3ZN9','KHP15','KHP102','KHP102M','KHP7025','KHP7026','KHP105','KHP109','C688','CU-Heet'}, ...
    [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]);

% 5) Neue numerische Spalte vorbereiten
numericLegierung = zeros(height(T),1);

% Fallback-Wert, falls unbekannte Legierung
fallbackValue = 1;
letzterWert = fallbackValue;

% 6) Schleife über Zeilen, Text -> Zahlenwert
for i = 1:height(T)
    % Aus dem Cell-Array mit {} auslesen
    currentValue = T.Legierung{i};

    % Überflüssige Anführungszeichen/Leerzeichen entfernen (sicherheitshalber)
    currentValue = strrep(currentValue, '"', '');
    currentValue = strtrim(currentValue);

    % Prüfen, ob im Mapping
    if isKey(legMap, currentValue)
        letzterWert = legMap(currentValue);
    else
        % Falls unbekannt -> bleibe beim letzten Wert
        % letzterWert = fallbackValue; 
        % (nur falls jedes Mal forcieren notwendig ist)
    end

    numericLegierung(i) = letzterWert;
end

% 7) Numerische Spalte in T ersetzen
T.Legierung = numericLegierung;

% 8) Unerwünschte Spalte löschen "ChannelInfoFieldText ..."
%    Der genaue Variablenname kann in T oder opts.VariableNames stehen
if ismember('ChannelInfoFieldText__10_1____PDA_expression__', T.Properties.VariableNames)
    T.ChannelInfoFieldText__10_1____PDA_expression__ = [];
end

% 9) Inf-Werte -> NaN
T{:,:}(isinf(T{:,:})) = NaN;

% 10) NaN-Werte durch "nächsten" Wert ersetzen
T = fillmissing(T,'next');

% Kontrolle:
head(T)


%% Spalten aus der Tabelle T in eigene Variablen übertragen
time = T.time;
t_band = T.Banddicke;
vzst_durchfluss = T.Verzinnungsstation_Durchflussmenge;
vzst_ruecklauf = T.Verzinnungsstation_R_cklauftemperatur;
vzst_kuehl = T.Verzinnungsstation_K_hlleistung;
vzstFU_durchfluss = T.VerzinnungsstationFU_Durchflussmenge;
vzstFU_ruecklauf = T.VerzinnungsstationFU_R_cklauftemperatur;
power = T.Ofenleistung;
temp_wanne_aussen = T.WannentemperaturAu_en;
temp_wanne_innen = T.WannentemperaturInnen;
soll_temp = T.SolltemperaturUnterZinnbadrolle;
b_band = T.Bandbreite;
soll_temp_Still_Rezept = T.LaufprogrammOnlineZinnbad_AnlagenstartwertRezept;
soll_temp_Still_Korrigiert = T.LaufprogrammOnlineZinnbad_AnlagenstartwertKorrigiert;
airkniveOS = T.AbblasungOSLuftdruck;
airkniveUS = T.AbblasungUSLuftdruck;
speed = T.Bandgeschwindigkeit;
ist_temp_uZR = T.TemperaturUnterZinnbadrolle;
Ofen_1_aktiv = T.MONVerzinnungseinheitVerzinnungsstationOfen1Aktiv;
Ofen_2_aktiv = T.MONVerzinnungseinheitVerzinnungsstationOfen2Aktiv;
legierung = T.Legierung;
l_band = T.Bandl_ngeAufAbhaspel;
bandlauf_on = T.Bandtransport;
% soll_temp_start_korrigiert
traverse_down = T.HubtraverseUnten;

% Erstelle Arrays mit je einem Zeit- und einem Wertvektor
data_time = [time, time];
data_banddicke = [time, t_band];
data_vzst_durchfluss = [time, vzst_durchfluss];
data_vzst_ruecklauf = [time, vzst_ruecklauf];
data_vzst_kuehl = [time, vzst_kuehl];
data_vzstFU_durchfluss = [time, vzstFU_durchfluss];
data_vzstFU_ruecklauf = [time, vzstFU_ruecklauf];
data_power = [time, power];
data_wanne_temp_aussen = [time, temp_wanne_aussen];
data_wanne_temp_innen = [time, temp_wanne_innen];
data_soll_temp = [time, soll_temp];
data_bandbreite = [time, b_band];
data_soll_temp_still_rezept = [time, soll_temp_Still_Rezept];
data_soll_temp_still_korrigiert = [time, soll_temp_Still_Korrigiert];
data_airkniveOS = [time, airkniveOS];
data_airkniveUS = [time, airkniveUS];
data_speed = [time, speed];
data_ist_temp_uZR = [time, ist_temp_uZR];
data_ofen_1 = [time, Ofen_1_aktiv];
data_ofen_2 = [time, Ofen_2_aktiv];
data_legierung = [time, legierung];
data_bandlaenge = [time, l_band];
data_bandlauf_on = [time, bandlauf_on];
data_traverse_down = [time, traverse_down];