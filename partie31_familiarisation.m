clear all;
close all;


%%%%%%%%%%% Partie 3.1 familiarisation %%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Définir les paramètres du code convolutif
treillis = poly2trellis(3, [5, 7]); % Treillis du code convolutif (longueur de contrainte = 3, générateurs [5, 7] en octal)
message = [0 1 0 1 1 0 1]; % Message binaire d'entrée
code = convenc(message, treillis); % Encodage convolutif

tb=5;
decoded = vitdec(code,treillis,tb,'trunc','hard');

% Nombre d'états dans le treillis (2^(longueur de contrainte - 1))
nombreEtats = treillis.numStates;

% Obtenir les transitions et sorties à partir du treillis
transitions = treillis.nextStates + 1; % Transitions entre états (+1 pour indices MATLAB)
sorties = treillis.outputs; % Sorties associées aux transitions

% Définir le nombre de bits de sortie (nombre de générateurs dans [5, 7])
nombreBitsSortie = length([5, 7]); % Nombre de bits de sortie par transition

% Tracer le treillis
figure;
hold on;
title('Treillis du code convolutif'); % Titre du graphique
xlabel('Temps'); ylabel('État'); % Légendes des axes

% Parcourir les états et les transitions pour dessiner les flèches
for etatActuel = 1:nombreEtats
    for bitEntree = 0:1
        etatSuivant = transitions(etatActuel, bitEntree + 1); % État suivant
        bitsSortie = dec2bin(sorties(etatActuel, bitEntree + 1), nombreBitsSortie); % Bits de sortie
        
        % Tracer une flèche entre les états
        plot([0 1], [etatActuel etatSuivant], '-o', 'LineWidth', 1.5);
        text(0.5, (etatActuel + etatSuivant)/2, bitsSortie, ...
            'HorizontalAlignment', 'center', 'BackgroundColor', 'white'); % Étiquette des bits de sortie
    end
end
&&&&&&&&&&&&&&&&&&&&&&& la figure affichée est cellule élémentaire du treillis &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&& d'apres ce que on a compris d'un cours mis en ligne, on vous met ci joint le lien du cours(slider 11), une cellule élémentaire 
&&&&&&&& représente toutes les transitions possible. Du coup plus facile pour l'algorithme de viterbi.
&&&&&&&&
&&&&&&&& Lien : https://ensat.ac.ma/Portail/wp-content/uploads/2020/03/Modul25_TI_cours3.pdf
&&&&&&&&
&&&&&&&&
&&&&&&&& 
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

%Ajuster l'affichage
yticks(1:nombreEtats);
yticklabels(arrayfun(@(x) ['État ' num2str(x-1)], 1:nombreEtats, 'UniformOutput', false)); % Étiquettes des états
grid on;
hold off;

% Afficher le message original
disp('Message binaire d''entrée :');
disp(message);

% Afficher le code convolutif encodé
disp('Message encodé (code convolutif) :');
disp(code);

% Afficher le message décodé après utilisation de l'algorithme de Viterbi
disp('Message décodé (après le décodage avec Viterbi) :');
disp(decoded);

% Afficher le nombre d'erreurs entre les données originales et les données décodées
disp('Nombre d''erreurs entre le message original et le message décodé :');
disp(biterr(message, decoded));
