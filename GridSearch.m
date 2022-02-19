r2O = r2; %R� original
alphaO = par(1); %Alpha Original
betaO = par(2); %Beta Original
passo = 0.1; %Passo inicial de 0.1

R2_grad = 1;
alphaI = alphaO;
betaI = betaO;
it=0;

while(norm(R2_grad) >= 0.0001) %Precis�o de R� definido aqui
    it=it+1;
    [R2, a, b] = procurar(passo, edges, alphaI, betaI, VPreal, theta); %criar a matriz R2 de R� em torno de alpha e beta iniciais, a � o vetor dos alphas e b � o vetor de betas para gerar surf(a, b, R2)
    r2I = max(max(R2)); %r2I � o maior valor da matriz R2 obtido iterativamente
    [i1, i2] = find(R2==r2I); %encontra a posi��o do maior R� de R2
    alphaI = a(i1); %encontra o valor de alpha para a posi��o do maior R�
    betaI = b(i2); %encontra o valor de beta para a posi��o do maior R�
    [fx, fy] = gradient(R2,passo); %encontra os gradientes de todos os pontos da matrix R2
    R2_grad = [fx(i1) fy(i2)]; %encontra a taxa de varia��o no ponto m�ximo de R2
    
    
    
    
    
    if(passo==0.1)
        f2 = figure;
        surf(a,b,R2') %Plot de superf�cie. A FORMA COMO R� � CRIADO � TRANSPOSTA, PREENCHENDO COLUNA POR COLUNA EM VEZ DE LINHA POR LINHA, por isso faz-se necess�rio a transposi��o para se obter os pontos desejados no plot
        hold on
        pI = scatter3(alphaI ,betaI, R2(i1, i2)+0.001, 'LineWidth', 4) %ponto obtido
        pO = scatter3(alphaO, betaO, r2O+0.001, 'LineWidth', 4) %ponto inicial
        title(['Grid Search - ' , place], 'FontSize', 14)
        legend([pI pO],{'M�ximo R�','Ponto Inicial'}, 'FontSize', 14, 'Location', 'northeast')
        %ALPHA=a;
        %BETA=b;
    end
    if(passo==0.1*0.1)
       figure(f2)
       hold on
       ALPHA=a;
       BETA=b;
       surf(ALPHA, BETA, R2', 'HandleVisibility','off')
    end
    
    
            
    xlabel('\beta', 'FontWeight', 'bold', 'FontSize', 14) %INVERTIDO ALPHA POR BETA POR CAUSA DO ARTIGO
    ylabel('\alpha', 'FontWeight','bold', 'FontSize', 14) %INVERTIDO ALPHA POR BETA POR CAUSA DO ARTIGO
    zlabel('R�', 'FontWeight', 'bold', 'FontSize', 14)
    title(colorbar,'R�', 'FontSize', 14)
    xticks(0:0.5:5)
    yticks(0:0.5:5)
    
    passo = passo*0.1; %incrementa passo para pr�xima itera��o se condi��o n�o estabelicida
        
end


function [y, a, b] = procurar(passo, edges, alpha, beta, VPreal, theta)
i1 = 0;
i2 = 0;

for ia = (alpha-(10*passo)):passo:(alpha+(10*passo))
    i1 = i1+1;
    for ib = (beta-(10*passo)):passo:(beta+(10*passo))
        i2 = i2+1;
        
        a1 = wblcdf(edges(find(edges>0)), ia, ib);
        a2 = wblcdf(edges(find(edges<15)), ia, ib);
        VPestimated = (a1-a2)*(1-theta);
        VPestimated(1)=theta+VPestimated(1);
        r = corrcoef(VPreal, VPestimated);
        R2(i1, i2) = r(2)^2;
    end
    i2 = 0;
end
i1 = 0;

y = R2;
a = (alpha-(10*passo)):passo:(alpha+(10*passo));
b = (beta-(10*passo)):passo:(beta+(10*passo));
end



