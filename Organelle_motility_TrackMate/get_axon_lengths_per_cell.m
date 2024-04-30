addpath('/Volumes/Emily_htt_2/AC_codes_epmodified_20200603/');

for k_choose = 4
    if k_choose == 1    % 30Q 
        cd('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_mats/');
        fl='*.mat';
        save_dir='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_dir/';
    elseif k_choose == 2   % 45Q 
        cd('/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_lyso_mats/');
        fl='*.mat';
        save_dir='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_lyso_dir/';
    elseif k_choose == 3    % 65Q 
        cd('/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_lyso_mats/');
        fl='*.mat';
        save_dir='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_lyso_dir/';
    elseif k_choose == 4    % 81Q 
        cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_mats/');
        fl='*.mat';
        save_dir='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_dir/';
    elseif k_choose == 5    % 30Q+IFg 
        cd('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_IFg_mats/');
        fl='*.mat';
        save_dir='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_IFg_dir/';
    elseif k_choose == 6    % 81Q+IFg
        cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_IFg_mats/');
        fl='*.mat';
        save_dir='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_IFg_dir/';
    end
    
    fls=dir(fl);
    axon_length=cell(numel(fls),1);
    filename=cell(numel(fls),1);

    for k=1:numel(fls)
        load(fls(k).name);
        axon_length{k}=ab.axon_length;
        filename{k}=fls(k).name;
    end
    save([save_dir, 'axon_length_per_cell'],'filename','axon_length');

end