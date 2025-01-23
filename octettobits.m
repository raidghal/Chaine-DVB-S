function bits = octettobits(octets)
    % Convertir les nombres décimaux en séquences de 8 bits
    % decimal_vector : vecteur de nombres décimaux (valeurs entre 0 et 255)    
    % Initialiser le vecteur de bits
    bits = zeros(1, 8 * length(octets));
    
    % Convertir chaque nombre décimal en 8 bits
    for i = 1:length(octets)
        decimal = octets(i);
        for j = 1:8
            bits((i-1)*8 + (9 - j)) = mod(decimal, 2); % Bit de poids faible en premier
            decimal = floor(decimal / 2); % Diviser par 2 pour passer au bit suivant
        end
    end
end