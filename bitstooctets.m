function octets = bitstooctets(bits)
    % Convertir une séquence de bits en nombres décimaux (octets)
    % bits : vecteur de bits (longueur doit être un multiple de 8)
    % decimal_vector : vecteur de nombres décimaux (valeurs entre 0 et 255)
    
    % Vérifier que la longueur des bits est un multiple de 8
    if mod(length(bits), 8) ~= 0
        error('La longueur de la séquence de bits doit être un multiple de 8.');
    end
    
    % Regrouper les bits par 8
    bits_reshape = reshape(bits, 8, []).'; % Chaque ligne = 8 bits
    
    % Convertir chaque groupe de 8 bits en un nombre décimal
    octets = zeros(1, size(bits_reshape, 1)); % Initialiser le vecteur de sortie
    for i = 1:size(bits_reshape, 1)
        somme = 0;
        for j = 1:8
            somme = somme + bits_reshape(i, j) * 2^(8 - j); % Conversion binaire -> décimal
        end
        octets(i) = somme;
    end
end