function ChangeAxisLabel(FigNum,Size,NewXLabel,NewTitle)

figure(FigNum)
for ii = 1:Size(1)
    for jj = 1:Size(2)
            if ii == 1 
                nn = ii + jj - 1;
            else
                nn = ii + jj - 1;
            end
        subplot(Size(1),Size(2),nn)
        xlabel(NewXLabel,'Interpreter','latex','FontSize',14)
    end
end

sgtitle(NewTitle,'Interpreter','latex')


end