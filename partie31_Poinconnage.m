clear all 
close all 

&&&&&&&&& Partie31 Poinconnage &&&&&&&&&&&&&&&&&&&&


&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



%% Constantes d'entrée
M = 4; % Ordre de modulation
l = log2(M); % Nombre de bits par symbole
alpha = 0.35; % Roll-off du filtre SRRC (Square Root Raised Cosine)           
span = 10; % Longueur du filtre SRRC (en nombre de périodes symboles). 
Ns = 5; % Facteur de suréchantillonnage
N = 204;
K = 188;
nbits = 188*204*2; % Nombre total de bits à transmettre. 


% Génération de bits et mapping QPSK
bits = randi([0,1], 1, nbits); % Bits aléatoires

%% Code convolutif et poinçonnement

% Définition du treillis pour le codage convolutif 
trellis = poly2trellis(7, [171 133]); 

% Matrice de poinçonnement pour réduire le taux du code R=1/2 -> R=3/4
puncmat = [1 1 0 1]; %baisse de la redondance
code = convenc(bits, trellis); % Codage convolutif des bits sans poinçonnement

% Mapping QPSK 
symboles = (1 - 2 * code(1:2:end)) + 1i * (1 - 2 * code(2:2:end)); 


% Suréchantillonnage
signal = kron(symboles, [1 zeros(1, Ns-1)]); % Ajout de zéros pour suréchantillonner



%% Filtre de mise en forme

% SPAN est la durée du filtre en nombre période symbole
% Les filtres générés par la fonction rcosdesign ont une durée de SPAN
% périodes symboles, soit SPAN * SPS +1 périodes d'échantillonage SPS = Ns
% facteur de sur-échantillonage

% Génération du filtre en cosinus surélevé avec roll-off alpha
h = rcosdesign(alpha, span, Ns, 'sqrt');
signal_filtre = filter(h, 1, [signal zeros(1, length(h)-1)]); % Filtrage émission

%% Code convolutif et poinçonnement

% Codage convolutif avec poinçonnement
codepunct = convenc(bits,trellis,puncmat);

% Mapping QPSK
symbolespunct = (1 - 2 * codepunct(1:2:end)) + 1i * (1 - 2 * codepunct(2:2:end)); % Symboles complexes
signalpunct = kron(symbolespunct, [1 zeros(1, Ns-1)]); % Ajout de zéros pour suréchantillonner
signal_filtre_punct = filter(h, 1, [signalpunct zeros(1, length(h)-1)]); % Filtrage émission


%% Canal AWGN
TEBsoft = [];
TEB_theo = [];
EbN0 = -4:4;
EbN0lin=10.^(EbN0/10);
TEB_theo = qfunc(sqrt(2 * EbN0lin));


for i = 1:length(EbN0)
    %AWGN
    
    % Calcul de la puissance moyenne du signal filtré
    Px = mean(abs(signal_filtre).^2); 
    
    % Calcul de la variance du bruit en fonction de Eb/N0
    sigma = (Px * Ns) / (2 * l * EbN0lin(i)); 
    
    % Génération de bruit gaussien complexe (réel et imaginaire indépendants)
    bruit_reel = sqrt(sigma) * randn(size(signal_filtre)); 
    bruit_imag = sqrt(sigma) * randn(size(signal_filtre)); 
    bruit = bruit_reel + 1i * bruit_imag; 
    
    % Ajout du bruit au signal pour obtenir le signal bruité
    signal_bruite = signal_filtre + bruit; 


    % Réception
    % Filtrage du signal bruité avec un filtre en cosinus surélevé (SRRC)
    signal_recu_bruite = filter(h, 1, signal_bruite); 
    
    % Extraction des symboles échantillonnés au rythme symbolique
    signal_echantillonne_bruite = signal_recu_bruite(length(h):Ns:end); 

    
    % Démapping QPSK
    bitsrecus_reels = real(signal_echantillonne_bruite); % Partie réelle
    bitsrecus_imag = imag(signal_echantillonne_bruite); % Partie imaginaire
    
    % Reconstruction des bits
    bitsrecus = zeros(1, length(code));
    bitsrecus(1:2:end) = bitsrecus_reels;
    bitsrecus(2:2:end) = bitsrecus_imag;
    
    tb = 30; % code 1/2 donc 5*(7-1) matlab
    
    %% Decodage soft sans poinconnement
    % sans poinconnement
    decodedsoft = vitdec(bitsrecus,trellis,tb,'trunc','unquant');

    % Calcul du TEB soft
    ecartsoft = sum(bits ~= decodedsoft); % Différence entre bits transmis et reçus
    TEBsoft(i) = ecartsoft/nbits;
    
    %% Decodage soft et hard avec poinconnement
    %AWGN poinçonnement
    Px = mean(abs(signal_filtre_punct).^2);
    sigma = (Px*Ns) / (2*l*EbN0lin(i));
    bruit_reel = sqrt(sigma) * randn(size(signal_filtre_punct));
    bruit_imag = sqrt(sigma) * randn(size(signal_filtre_punct));
    bruit = bruit_reel + 1i *bruit_imag;
    signal_bruite_punct = signal_filtre_punct + bruit;
    % Réception
    signal_recu_bruite_punct = filter(h, 1,signal_bruite_punct);
    signal_echantillonne_bruite_punct = signal_recu_bruite_punct(length(h):Ns:end); % Extraction des symboles
     % Démapping QPSK
    bitsrecus_reels_punct = real(signal_echantillonne_bruite_punct); % Partie réelle
    bitsrecus_imag_punct = imag(signal_echantillonne_bruite_punct); % Partie imaginaire
    
    % Reconstruction des bits
    bitsrecus_punct = zeros(1, length(codepunct));
    bitsrecus_punct(1:2:end) = bitsrecus_reels_punct;
    bitsrecus_punct(2:2:end) = bitsrecus_imag_punct;
    
    tb = 30; % code 1/2 donc 5*(7-1) matlab
    
    %decodage poinconnement soft 
    decodedsoft_p = vitdec(bitsrecus_punct,trellis,tb,'trunc','unquant',puncmat);
    
    ecartsoft_p = sum(bits ~= decodedsoft_p); % Différence entre bits transmis et reçus
    TEBsoftp(i) = ecartsoft_p/nbits;

   
    
        
end


%% Affichage des résultats
figure;

% Tracé des courbes
semilogy(EbN0, TEBsoft, 'b-', 'LineWidth', 2); % Courbe soft

hold on;
semilogy(EbN0, TEBsoftp, 'g-', 'LineWidth', 2); % Courbe soft poinconne

% Personnalisation de l'affichage
xlabel('SNR [dB]');
ylabel('TEB');
legend( 'Soft', 'Soft avec poinçonnement', 'FontSize', 10, 'Location', 'southwest');
title('Tracé des TEBs avec et sans poinçonnement en fonction de Eb/N0');
% Ajustement des limites de l'axe des Y pour inclure une précision jusqu'à 10^-3
ylim([1e-3, 1]); % Axe Y entre 10^-3 et 1
grid on;

hold off;



