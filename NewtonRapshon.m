tic
%close all
%clear all
%clc

load('serie_dados') %Carrega os dados da Raiz
%Descomentar as 3 seguintes linhas se os dados forem de CuritibaUTFPR
% CuritibaUTFPR = load('CuritibaUTFPR', 'vveloc'); 
% CuritibaUTFPR = CuritibaUTFPR.('vveloc');
% CuritibaUTFPR = CuritibaUTFPR-0.4;
edges = -0.05:0.1:15.05;

%Substituir o nome dos 3 próximos vetores pelos nomes dos arquivos na raiz
par = wblfit(serie_dados(find(serie_dados>=0.1))); %Valores iniciais para alpha e beta estimados pelo método da máxima verossilmilhança
Hist = histogram(serie_dados, edges, 'Normalization', 'pdf'); %Hist carrega várias informações como os proprios dados (Hist.Data), os valores usados no histograma (Hist.Values) e o total de valores medidos para cada bin (Hist.BinCounts)
axis([-0.05 8 0 0.6])
xticks(0:1:15)
xlabel('Wind speed (m/s)', 'FontSize', 14)
ylabel('Frequency', 'FontSize', 14)
title(sprintf('Curitiba - UTFPR'), 'FontSize', 14)
grid
hold on

P = 0.00001; %Intervalo entra um alpha[n] e alpha[n-1] tal como para os betas
epsilon = 0.000001; %Critério de parada a ser comparado com as equações A=d(R²)/dalpha e B=(R²)/dbeta
it = 1;
format long %mostrar mais casas decimais
[R2, theta] = rsquared(Hist, par) %função rsquared calcula R² em função do histograma e dos parâmetros alpha e beta
pdf = wblpdf(0:0.001:15, par(1), par(2))*(1-theta);
plot(0:0.001:15, pdf, 'LineWidth', 3, 'Color', [200/255 90/255 0/255])
St = stem(0, Hist.Values(1), 'LineWidth', 3, 'Color', [255/255 0/255 0/255], 'HandleVisibility','off', 'Marker', '^')
txtr2O = ['R² = ', num2str(R2)];
text(2.5, 0.35, txtr2O, 'Color', [200/255 90/255 0/255], 'FontWeight', 'bold', 'FontSize', 18)

while(1)
    
    alpha = [par(1)-P par(1) par(1)+P]; % cria dois pontos em torno do ponto inicial de alpha com o passo P
    beta = [par(2)-P par(2) par(2)+P]; % cria dois pontos em torno do ponto inicial de beta com o passo P
    %% Alpha variando e Beta constante
    
    A = ( rsquared(Hist, [alpha(2) beta(2)]) - rsquared(Hist, [alpha(1) beta(2)]) )/P; %A=d(R²)/dalpha por Backward Euler
    Ait(it) = A; %cria vetor de A em função da iteração
    %% Beta variando e Alpha constante
        
    B = ( rsquared(Hist, [alpha(2) beta(2)]) - rsquared(Hist, [alpha(2) beta(1)]) )/P; %A=d(R²)/dbeta por Backward Euler
    Bit(it) = B; %cria vetor de B em função da iteração
    %% Teste de convergência
    if(epsilon > max(abs(A), abs(B)))
        break;
    end
    
    %% Calcula os R² em torno do ponto incial exceto os pontos que não são utilizados
    MR2 = [rsquared(Hist, [alpha(1) beta(1)]) rsquared(Hist, [alpha(1) beta(2)])                    0
           rsquared(Hist, [alpha(2) beta(1)]) rsquared(Hist, [alpha(2) beta(2)]) rsquared(Hist, [alpha(2) beta(3)])
                         0                    rsquared(Hist, [alpha(3) beta(2)])                    0              ];
    
    %% Matriz Jacobiana
    H = (1/(P^2))*( MR2(3,2) - 2*MR2(2,2) + MR2(1,2) ); %dA/d\alpha por Second-order central
    N = (1/(P^2))*( MR2(2,2) - MR2(2,1) - MR2(1,2) + MR2(1,1) ); %dA/d\beta por Backward Euler 2 vezes
    M = (1/(P^2))*( MR2(2,2) - MR2(1,2) - MR2(2,1) + MR2(1,1) ); %dB/d\alpha por Backward Euler 2 vezes
    L = (1/(P^2))*( MR2(2,3) - 2*MR2(2,2) + MR2(2,1) ); %dB/d\beta por Second-order central
    J = [H N; M L];
    
    %% Atualização de parâmetros x[n+1] = x[n] - inv(J)*g(x)
    gx = [A; B];
    par = (par'-inv(J)*gx)';
    alphait(it) = par(1);
    betait(it) = par(2);
    it = it+1;
%     R2 = rsquared(Hist, par) %calcula R² para cada iteração e mostra no Command Window
end
R2 = rsquared(Hist, par)
pdf = wblpdf(0:0.001:15, par(1), par(2))*(1-theta);
plot(0:0.001:15, pdf, 'LineWidth', 3, 'Color', [110/255 0/255 190/255])
stem(0, Hist.Values(1), 'LineWidth', 3, 'Color', [110/255 0/255 190/255], 'HandleVisibility','off')
txtr2O = ['R² = ', num2str(R2)];
text(2.5, 0.3, txtr2O, 'Color', [110/255 0/255 190/255], 'FontWeight', 'bold', 'FontSize', 18)
toc
legend({'Histogram', 'Maximum Likelihood', 'Newton Raphson'},'FontSize', 12)
set(gcf, 'Position', [100 100 700 500])
ffpdf = fullfile('Imagens', 'CuritibaUTFPF_PDF.eps');
saveas(gcf, ffpdf)

% figure
% plot(CuritibaUTFPR)
% ylabel('Wind speed (m/s)', 'FontSize', 14)
% xlabel('Acquisition', 'FontSize', 14)
% title(sprintf('Curitiba - UTFPR'), 'FontSize', 14)
% axis([0 7.5e4 0 9])
% set(gcf, 'Position', [100 100 700 500])
% ffpdf = fullfile('Imagens', 'CuritibaUTFPF_Time.eps');
% saveas(gcf, ffpdf)
% ffpdf = fullfile('Imagens', 'CuritibaUTFPF_Time.png');
% saveas(gcf, ffpdf)

% figure
% plot(CuritibaUTFPR)
% ylabel('Wind speed (m/s)', 'FontSize', 14)
% xlabel('Acquisition', 'FontSize', 14)
% title(sprintf('Curitiba - UTFPR'), 'FontSize', 14)
% axis([0 1000 0 9])
% set(gcf, 'Position', [100 100 700 500])
% ffpdf = fullfile('Imagens', 'CuritibaUTFPF_Time2.eps');
% saveas(gcf, ffpdf)
% ffpdf = fullfile('Imagens', 'CuritibaUTFPF_Time2.png');
% saveas(gcf, ffpdf)

% figure
% plot(1:it, Ait)
% xlabel('Iterações')
% ylabel('$\frac{dR^2(\alpha,\beta)}{d\alpha}$', 'Interpreter', 'Latex', 'FontSize', 20)
% figure
% plot(1:it, Bit)
% xlabel('Iterações')
% ylabel('$\frac{dR^2(\alpha,\beta)}{d\beta}$', 'Interpreter', 'Latex', 'FontSize', 20)
% figure
% plot(1:it-1, alphait)
% xlabel('Iterações')
% ylabel('$\alpha$', 'Interpreter', 'Latex', 'FontSize', 20)
% figure
% plot(1:it-1, betait)
% xlabel('Iterações')
% ylabel('$\beta$', 'Interpreter', 'Latex', 'FontSize', 20)
