clear all;
close all;


&&&&&&&&&&&& Partie 31 Hard et Soft &&&&&&&&&&&&&&&&&&&&&&&&&&&&&


&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&




%% Paramètres d'entrée
M = 4; % Ordre de modulation (QPSK)
l = log2(M); % Nombre de bits par symbole
nbits = 188 * 8; % Nombre total de bits à transmettre
alpha = 0.35; % Roll-off du filtre SRRC
span = 10; % Longueur du filtre SRRC (en nombre de périodes symboles)
Ns = 5; % Facteur de suréchantillonnage
tb = 30; % Taille de la trace-back (code 1/2 donc 5*(7-1))

%% Génération des bits et codage convolutif
bits = randi([0, 1], 1, nbits); % Génération des bits aléatoires

trellis = poly2trellis(7, [171 133]); % Génération du trellis

code = convenc(bits, trellis); % Codage convolutif

% Mapping QPSK
symboles = (1 - 2 * code(1:2:end)) + 1i * (1 - 2 * code(2:2:end)); % Symboles complexes

% Suréchantillonnage

signal = kron(symboles, [1 zeros(1, Ns - 1)]);

% Filtre de mise en forme
% Racine de cosinus surélevé alpha = 0.35
h = rcosdesign(alpha, span, Ns, 'sqrt');

% Signal en sortie 
signal_filtre = filter(h, 1, [signal zeros(1, length(h) - 1)]);

%% Canal AWGN
% Définition du rapport signal sur bruit par bit (Eb/N0) en dB
EbN0 = -4:4; % Rapport Eb/N0 en dB (de -4 dB à 4 dB)
EbN0lin = 10.^(EbN0 / 10); % Conversion de Eb/N0 de l'échelle logarithmique à l'échelle linéaire



% Initialisation des vecteurs pour stocker les TEB obtenus par décodage hard et soft
TEBhard = [];
TEBsoft = [];


for i = 1:length(EbN0) 
    
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

    %% Décodage hard 
    % Les bits sont récupérés en comparant les parties réelle et imaginaire des symboles à 0
    bitsrecus_reels = real(signal_echantillonne_bruite) < 0; % Partie réelle
    bitsrecus_imag = imag(signal_echantillonne_bruite) < 0; % Partie imaginaire
    
    % Reconstruction des bits décodés à partir des parties réelle et imaginaire
    bitsrecus = zeros(1, length(code)); 
    bitsrecus(1:2:end) = bitsrecus_reels; % Bits récupérés à partir des parties réelles
    bitsrecus(2:2:end) = bitsrecus_imag; % Bits récupérés à partir des parties imaginaires
    
    % Décodage convolutif de viterbi hard
    decodedhard = vitdec(bitsrecus, trellis, tb, 'trunc', 'hard'); 

    %% Décodage soft 
    % Les parties réelle et imaginaire des symboles sont directement utilisées
    bitsrecus_reels_soft = real(signal_echantillonne_bruite); 
    bitsrecus_imag_soft = imag(signal_echantillonne_bruite); 
    
    % Reconstruction des bits décodés en mode soft
    bitsrecus_soft = zeros(1, length(code)); 
    bitsrecus_soft(1:2:end) = bitsrecus_reels_soft; % Partie réelle utilisée sans quantification
    bitsrecus_soft(2:2:end) = bitsrecus_imag_soft; % Partie imaginaire utilisée sans quantification
    
    % Décodage convolutif de viterbi soft 
    decodedsoft = vitdec(bitsrecus_soft, trellis, tb, 'trunc', 'unquant'); 

    % Calcul des Taux d'Erreur Binaire (TEB)
    % Comparaison entre les bits transmis et les bits décodés pour les modes hard et soft
    TEBhard(i) = sum(bits ~= decodedhard) / nbits; % TEB pour le décodage hard
    TEBsoft(i) = sum(bits ~= decodedsoft) / nbits; % TEB pour le décodage soft
end


%% Affichage des résultats
figure;

% Tracé des TEB hard et soft
semilogy(EbN0, TEBhard, 'g-', 'LineWidth', 2); % Courbe hard
hold on;
semilogy(EbN0, TEBsoft, 'b-', 'LineWidth', 2); % Courbe soft

xlabel('SNR [dB]', 'FontSize', 12);
ylabel('TEB', 'FontSize', 12);

legend('TEB Hard', 'TEB Soft', 'FontSize', 10, 'Location', 'southwest');

title('Tracé des TEBs :  Hard et Soft Decoding', 'FontSize', 14);

% Ajustement des limites de l'axe des Y pour inclure une précision jusqu'à 10^-3
ylim([1e-3, 1]); % Axe Y entre 10^-3 et 1
grid on;

hold off;

