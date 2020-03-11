function oqs_save_fig(fig, fn)

savefig(sprintf('%s.fig', fn));
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', sprintf('%s.pdf', fn));

end
