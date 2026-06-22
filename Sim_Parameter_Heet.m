    % Zinnbad (reines Zinn flüssig, ab 231,93°C)
        rohSn = 6960; % kg/m³ Dichte bei 270°C
        cpSn = 230; % J/(kg K) spezifische Wärmekapazität(flüssig) ab 231,11°C

    % Störgrößen
        % Umwelt
        TU = 35; % °C Umgebungstemperatur    
 
    % Zinnbad
        lB = 1.44; % m Länge ZinnBad
        bB = 0.89; % m Breite ZinnBad
        hB = 0.86; % m Höhe ZinnBad
        VB = lB*bB*hB; % m³, Volumen
        mB = VB*rohSn; % kg, Masse (Zinn)
        CB = mB*cpSn; % J/K, Effektive Leistungsaufnahme gesamtes Bad 