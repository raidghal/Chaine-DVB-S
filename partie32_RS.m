clear all;
close all;

%%%%%%%%%%%% Partie32 Reed Solomon %%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Constantes d'entrée
M = 4; % Ordre de modulation
l = log2(M); % Nombre de bits par symbole
alpha = 0.35; % Roll-off du filtre SRRC
span = 10; % Longueur du filtre SRRC (en périodes symboles)
Ns = 5; % Facteur de suréchantillonnage
N = 204;
K = 188;
nbits = 188 *N*2; % Nombre total de bits à transmettre





%% Génération de bits et codage convolutif avec et sans RS
bits = randi([0, 1], 1, nbits); % Bits aléatoires

% Création du treillis 
trellis = poly2trellis(7, [171 133]);

% Matrice de poinconnage 
puncmat = [1 1 0 1];

%% Avec RS

% Mapping des bits en symboles QPSK


% Code externe Reed Solomon
encoder_rs = comm.RSEncoder(N, K, BitInput=true);
% Ajout du codage externe de Reed-Solomon
bits_rs = step(encoder_rs, bits.').';
code_rs_conv = convenc(bits_rs, trellis, puncmat);

% Mapping QPSK Reed solomon

symboles_rs = (1 - 2 * code_rs_conv(1:2:end)) + 1i * (1 - 2 * code_rs_conv(2:2:end));
% Suréchantillonnage
signal_rs = kron(symboles_rs, [1 zeros(1, Ns - 1)]);
h = rcosdesign(alpha, span, Ns, 'sqrt');
signal_rs_filtre = filter(h, 1, [signal_rs zeros(1, length(h) - 1)]);

%% Sans RS
code_conv = convenc(bits, trellis, puncmat);
% Mapping QPSK
symboles = (1 - 2 * code_conv(1:2:end)) + 1i * (1 - 2 * code_conv(2:2:end));
% Suréchantillonnage
signal = kron(symboles, [1 zeros(1, Ns - 1)]);
% Filtre de mise en forme
h = rcosdesign(alpha, span, Ns, 'sqrt');
signal_filtre = filter(h, 1, [signal zeros(1, length(h) - 1)]);

%% Canal AWGN
EbN0 = -4:4; % Rapport Eb/N0 en dB
EbN0lin = 10.^(EbN0 / 10);
TEB_theo = qfunc(sqrt(2 * EbN0lin)); % TEB théorique pour la QPSK
TEBsoft_rs = [];
TEBsoft = [];

for i = 1:length(EbN0)
    %% Avec REed solomon
    % Bruit AWGN RS
    Px_rs = mean(abs(signal_rs_filtre).^2);
    sigma_rs = (Px_rs * Ns) / (2 * l * EbN0lin(i));
    bruit_reel_rs = sqrt(sigma_rs) * randn(size(signal_rs_filtre));
    bruit_imag_rs = sqrt(sigma_rs) * randn(size(signal_rs_filtre));
    bruit_rs = bruit_reel_rs + 1i * bruit_imag_rs;
    
    % Signal bruité
    signal_rs_bruite = signal_rs_filtre + bruit_rs;
        
    % Réception et échantillonnage

    signal_rs_recu = filter(h, 1, signal_rs_bruite);
    signal_rs_echantillonne = signal_rs_recu(length(h):Ns:end);
    
    bitsrecus_rs_reel = real(signal_rs_echantillonne);
    bitsrecus_rs_imag = imag(signal_rs_echantillonne);

    % Reconstruction des bits
    bitsrecus_rs = zeros(1, length(code_rs_conv));
    bitsrecus_rs(1:2:end) = bitsrecus_rs_reel;
    bitsrecus_rs(2:2:end) = bitsrecus_rs_imag;

    tb = 30; % Longueur de la fenêtre de traceback
    decodedsoft_rs = vitdec(bitsrecus_rs, trellis, tb, 'trunc', 'unquant', puncmat);
    
    % Decodage Reed-Salomon RS(204,188)
    decoder_rs = comm.RSDecoder(N, K,BitInput=true);
    
    decodedsoft_rs_final = step(decoder_rs, decodedsoft_rs.');
    decodedsoft_rs_final = decodedsoft_rs_final.';

    % Calcul des TEBs
    ecartsoft_rs = sum(bits ~= decodedsoft_rs_final);
    TEBsoft_rs(i) = ecartsoft_rs / nbits;
    %% Sans Reeed solomon
    % Bruit AWGN
    Px = mean(abs(signal_filtre).^2);
    sigma = (Px * Ns) / (2 * l * EbN0lin(i));
    bruit_reel = sqrt(sigma) * randn(size(signal_filtre));
    bruit_imag = sqrt(sigma) * randn(size(signal_filtre));
    bruit = bruit_reel + 1i * bruit_imag;
    
    % Signal bruité
    signal_bruite = signal_filtre + bruit;

    % Réception et échantillonnage
    signal_recu = filter(h, 1, signal_bruite);
    signal_echantillonne = signal_recu(length(h):Ns:end);

    bitsrecus_reel = real(signal_echantillonne) ;
    bitsrecus_imag = imag(signal_echantillonne) ;

    % Reconstruction des bits

    bitsrecus = zeros(1, length(code_conv));
    bitsrecus(1:2:end) = bitsrecus_reel;
    bitsrecus(2:2:end) = bitsrecus_imag;

    % Décodage convolutif
    
    decodedsoft = vitdec(bitsrecus, trellis, tb, 'trunc', 'unquant', puncmat);

    

    ecartsoft = sum(bits ~= decodedsoft);
    TEBsoft(i) = ecartsoft / nbits;
end

%% Affichage des résultats
figure;
semilogy(EbN0, TEBsoft_rs, 'g-', 'LineWidth', 2);
hold on;

semilogy(EbN0, TEBsoft, 'b-', 'LineWidth', 2);
xlabel('SNR [dB]', 'FontSize', 12);
ylabel('TEB', 'FontSize', 12);
legend( 'Soft RS', 'Soft Sans RS', 'FontSize', 10, 'Location', 'southwest');
title('Tracé des TEBs Soft RS et Soft Sans RS', 'FontSize', 14);
ylim([1e-3, 1]); % Axe Y entre 10^-3 et 1
grid on;

hold off;
