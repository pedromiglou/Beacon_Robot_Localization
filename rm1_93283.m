function rm1_93283(N,Dt,r,L,Vn,Wn,v)
    % N - número de faróis (4)
    % Dt - intervalo de tempo de amostragem dos sensores (1 s)
    % r - raio das rodas dos robôs (0.15 m)
    % L - separação/afastamento das rodas conforme o modelo cinemático (1 m)
    % Vn - incerteza (sigma) na velocidade linear a impor ao robô (input) (0.1 m/s)
    % Wn - incerteza (sigma) na velocidade angular a impor ao robô (input) (0.1 m/s)
    % v - velocidade linear usada para calcular a segmentação da trajetória (5 m/s)

    if nargin == 0
        N = 4;
        Dt = 1;
        r = 0.15;
        L = 1;
        Vn = 0.1;
        Wn = 0.1;
        v= 5;
    end

    P = [0 0 0];

    B = BeaconDetection(N, P);
    
    YQ = pchip([0 B.X], [0 B(1).Y B(2).Y B(3).Y B(4).Y], 0:(v*Dt):B(4).X);

    plot(0:(v*Dt):B(4).X, YQ);

    X = 0:(v*Dt):B(4).X;
    Y = YQ;
    thetas = [0];

    for i=2:length(X)-1
        theta = atan2(Y(i+1)-Y(i), X(i+1)-X(i));
        thetas = [thetas theta];
    end

    thetas= [thetas thetas(end)];

    X
    Y
    thetas
end
