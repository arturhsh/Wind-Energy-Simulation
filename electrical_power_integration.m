electrical_power = out.electrical_power;
mechanical_power = out.mechanical_power;

electrical_energy = trapz(out.electrical_power);
mechanical_energy = trapz(out.mechanical_power);
opt_energy = trapz(out.opt_power);

ne = electrical_energy/opt_energy
nm = mechanical_energy/opt_energy

%%Fazer varredura entre 0.05 e 0.2 (4 valores: 0, 07, 14, 21)
%%Momento de inercia: Pegar metade e dobro
%%Revisar artigo procurando mais detalhes
%%Kopt: pegar 90% e 110% do valor
%%Adicionar ruído branco no sinal do vento que vai para o gerador (depois
%%do integrador do modelo) 
