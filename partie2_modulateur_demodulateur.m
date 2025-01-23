clear all;
close all;


%%%%%%%%%%%%%%%%% Partie 2 Implantation du modulateur/demodulateur %%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constantes d'entrée
M = 4; % Ordre de modulation
l = log2(M); % Nombre de bits par symbole
nbits = 188*8; % Nombre total de bits à transmettre. 
alpha = 0.35; % Roll-off du filtre SRRC (Square Root Raised Cosine)           
span = 10; % Longueur du filtre SRRC (en nombre de périodes symboles). 
Ns = 5; % Facteur de suréchantillonnage
        
% Génération de bits et mapping QPSK
bits = randi([0,1], 1, nbits); % Bits aléatoires
symboles = (1 - 2 * bits(1:2:end)) + 1i * (1 - 2 * bits(2:2:end)); % Symboles complexes

% Suréchantillonnage
signal = kron(symboles, [1 zeros(1, Ns-1)]); % Ajout de zéros pour suréchantillonner
% Filtre de mise en forme
% SPAN est la durée du filtre en nombre période symbole
% Les filtres générés par la fonction rcosdesign ont une durée de SPAN
% périodes symboles, soit SPAN * SPS +1 périodes d'échantillonage SPS = Ns
% facteur de sur-échantillonage
h = rcosdesign(alpha, span, Ns, 'sqrt');
signal_filtre = filter(h, 1, [signal zeros(1, length(h)-1)]); % Filtrage émission

%Canal AWGN
TEB = [];
TEB_theo = [];
EbN0 = -4:4;
EbN0lin=10.^(EbN0/10);
for i = 1:length(EbN0)
    %AWGN
    Px = mean(abs(signal_filtre).^2);
    sigma = (Px*Ns) / (2*l*EbN0lin(i));
    bruit_reel = sqrt(sigma) * randn(size(signal_filtre));
    bruit_imag = sqrt(sigma) * randn(size(signal_filtre));
    bruit = bruit_reel + 1i *bruit_imag;
    
    signal_bruite = signal_filtre + bruit;
    % Réception
    signal_recu_bruite = filter(h, 1,signal_bruite);
    signal_echantillonne_bruite = signal_recu_bruite(length(h):Ns:end); % Extraction des symboles
    
 
    % Démapping QPSK
    bitsrecus_reels = real(signal_echantillonne_bruite) < 0; % Partie réelle
    bitsrecus_imag = imag(signal_echantillonne_bruite) < 0; % Partie imaginaire
    
    % Reconstruction des bits
    bitsrecus = zeros(1, nbits);
    bitsrecus(1:2:end) = bitsrecus_reels;
    bitsrecus(2:2:end) = bitsrecus_imag;
    
    % Calcul du TEB
    ecart = sum(bits ~= bitsrecus); % Différence entre bits transmis et reçus
    TEB(i) = ecart/nbits;
end
TEB_theo = qfunc(sqrt(2 * EbN0lin));
figure(1);
semilogy(EbN0, TEB_theo,'r-')  
hold on
semilogy(EbN0, TEB,'b-')                                 
xlabel('SNR[dB]')                                    
ylabel('TEB');                                         
legend('Theorique', 'Simulé');
title('Tracé des TEBs en fonction de Eb/N0 ');
grid on;
hold off;
