clear all;
close all;

%% Constantes d'entrée
Bw = 40;
M = 4; % Ordre de modulation
l = log2(M); % Nombre de bits par symbole
Ts = 1; % Durée d’un symbole
Rs = 31.2; % Débit symbole (sym/s)
Rb = Rs * l; % Débit binaire en bits par seconde
alpha = 0.35; % Roll-off du filtre SRRC
Fe = 24000; % Fréquence d'échantillonnage en Hz (Fe >= (1+alpha)*Rs)
span = 10; % Longueur du filtre SRRC (en périodes symboles)
Ns = 5; % Facteur de suréchantillonnage
N = 204;
K = 188;
nbits = 188 *N*2; % Nombre total de bits à transmettre

%% Génération de bits
bits = randi([0, 1], 1, nbits); % Bits aléatoires

%% Codage RS
% Initialisation des encodeurs Reed-Solomon
encoder_rs = comm.RSEncoder(N, K, BitInput=true);

% Codage RS avec entrelacement
bits_rs = step(encoder_rs, bits.').'; % Codage RS
bits_rs_entrelace = convintrlv(bits_rs, 8, 3); % Entrelacement des bits RS

% Codage RS sans entrelacement
bits_rs_sans = step(encoder_rs, bits.').'; % Codage RS sans entrelacement



%% Codage convolutif
% Définition du treillis et poinçonnement
trellis = poly2trellis(7, [171 133]);
puncmat = [1 1 0 1];

% Codage convolutif avec entrelacement
code_rs_conv_entrelace = convenc(bits_rs_entrelace, trellis, puncmat);

% Codage convolutif sans entrelacement
code_rs_conv_sans = convenc(bits_rs_sans, trellis, puncmat);

%% Mapping QPSK
% Avec entrelacement
symboles_rs = (1 - 2 * code_rs_conv_entrelace(1:2:end)) + 1i * (1 - 2 * code_rs_conv_entrelace(2:2:end));

% Sans entrelacement
symboles_sans = (1 - 2 * code_rs_conv_sans(1:2:end)) + 1i * (1 - 2 * code_rs_conv_sans(2:2:end));

% Suréchantillonnage
signal_rs = kron(symboles_rs, [1 zeros(1, Ns - 1)]); % Avec entrelacement
signal_sans = kron(symboles_sans, [1 zeros(1, Ns - 1)]); % Sans entrelacement

%% Filtrage avec SRRC
h = rcosdesign(alpha, span, Ns, 'sqrt');
signal_rs_filtre = filter(h, 1, [signal_rs zeros(1, length(h) - 1)]); % Avec entrelacement
signal_sans_filtre = filter(h, 1, [signal_sans zeros(1, length(h) - 1)]); % Sans entrelacement

%% Canal AWGN
EbN0 = -4:4; % Rapport Eb/N0 en dB

% Initialisation des TEB
TEBsoft_rs = [];
TEBsoft_sans = [];

for i = 1:length(EbN0)
    %% Avec entrelacement
    % Bruit AWGN
    Px_rs = var(signal_rs_filtre);
    EbN0lin = 10^(EbN0(i) / 10); % Conversion en échelle linéaire
    sigma_rs = (Px_rs * Ns) / (2 * l * EbN0lin);
    bruit_reel_rs = sqrt(sigma_rs) * randn(size(signal_rs_filtre));
    bruit_imag_rs = sqrt(sigma_rs) * randn(size(signal_rs_filtre));
    bruit_rs = bruit_reel_rs + 1i * bruit_imag_rs;
    
    % Signal bruité
    signal_rs_bruite = signal_rs_filtre + bruit_rs;

    % Réception et échantillonnage
    signal_rs_recu = filter(h, 1, signal_rs_bruite);
    signal_rs_echantillonne = signal_rs_recu(length(h):Ns:end);

    % Décodage
    bitsrecus_rs_reel = real(signal_rs_echantillonne);
    bitsrecus_rs_imag = imag(signal_rs_echantillonne);
    bitsrecus_rs = zeros(1, length(code_rs_conv_entrelace));
    bitsrecus_rs(1:2:end) = bitsrecus_rs_reel;
    bitsrecus_rs(2:2:end) = bitsrecus_rs_imag;

    % Décodage convolutif
    tb = 30;
    decodedsoft_rs = vitdec(bitsrecus_rs, trellis, tb, 'trunc', 'unquant', puncmat);

    % Désentrelacement
    bits_desentrelacement = convdeintrlv(decodedsoft_rs, 8, 3);
    decodedsoft_rs_desentrelace = circshift(bits_desentrelacement, -8 * 3 * 7);

    % Décodage RS
    decoder_rs = comm.RSDecoder(N, K, BitInput = true);
    decodedsoft_rs_final = step(decoder_rs, decodedsoft_rs_desentrelace.');
    decodedsoft_rs_final = decodedsoft_rs_final.';

    % Calcul du TEB
    ecartsoft_rs = sum(bits ~= decodedsoft_rs_final);
    TEBsoft_rs(i) = ecartsoft_rs / nbits;

    %% Sans entrelacement
    % Bruit AWGN
    Px_sans = var(signal_sans_filtre);
    sigma_sans = (Px_sans * Ns) / (2 * l * EbN0lin);
    bruit_reel_sans = sqrt(sigma_sans) * randn(size(signal_sans_filtre));
    bruit_imag_sans = sqrt(sigma_sans) * randn(size(signal_sans_filtre));
    bruit_sans = bruit_reel_sans + 1i * bruit_imag_sans;

    % Signal bruité
    signal_sans_bruite = signal_sans_filtre + bruit_sans;

    % Réception et échantillonnage
    signal_sans_recu = filter(h, 1, signal_sans_bruite);
    signal_sans_echantillonne = signal_sans_recu(length(h):Ns:end);

    % Décodage
    bitsrecus_sans_reel = real(signal_sans_echantillonne);
    bitsrecus_sans_imag = imag(signal_sans_echantillonne);
    bitsrecus_sans = zeros(1, length(code_rs_conv_sans));
    bitsrecus_sans(1:2:end) = bitsrecus_sans_reel;
    bitsrecus_sans(2:2:end) = bitsrecus_sans_imag;

    % Décodage convolutif
    decodedsoft_sans = vitdec(bitsrecus_sans, trellis, tb, 'trunc', 'unquant', puncmat);

    % Décodage RS
    decodedsoft_sans_final = step(decoder_rs, decodedsoft_sans.');
    decodedsoft_sans_final = decodedsoft_sans_final.';

    % Calcul du TEB
    ecartsoft_sans = sum(bits ~= decodedsoft_sans_final);
    TEBsoft_sans(i) = ecartsoft_sans / nbits;
end
TEB_theo = qfunc(sqrt(2 * EbN0lin)); % TEB théorique pour la QPSK

%% Affichage des résultats
figure;
semilogy(EbN0, TEBsoft_rs, 'g', 'LineWidth', 2); % TEB avec entrelacement
hold on;
semilogy(EbN0, TEBsoft_sans, 'b', 'LineWidth', 2); % TEB sans entrelacement

xlabel('SNR [dB]', 'FontSize', 12);
ylabel('TEB', 'FontSize', 12);
legend('Avec Entrelacement', 'Sans Entrelacement', 'FontSize', 10, 'Location', 'southwest');
title('Comparaison des TEB avec et sans entrelacement', 'FontSize', 14);
ylim([1e-3, 1]); % Limite de l'axe Y
grid on;
hold off;
