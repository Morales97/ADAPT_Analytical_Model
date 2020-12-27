%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT BER to Eb/N0 relation for different MCS
% Find Eb/N0 that gives 1e-6 BER
% Daniel Morales - May 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EbNodB=0:1:20;
BER_BPSK = berawgn(EbNodB, 'psk', 2, 'nondiff');
BER_QPSK = berawgn(EbNodB, 'psk', 4, 'nondiff');
BER_8PSK = berawgn(EbNodB, 'psk', 8, 'nondiff');
BER_16QAM = berawgn(EbNodB, 'qam', 16);
BER_64QAM = berawgn(EbNodB, 'qam', 64);


semilogy(EbNodB,BER_BPSK, 'DisplayName', 'BPSK')
hold on
semilogy(EbNodB,BER_QPSK, 'DisplayName', 'QPSK')
semilogy(EbNodB,BER_8PSK, 'DisplayName', '8-PSK')
semilogy(EbNodB,BER_16QAM, 'DisplayName', '16-QAM')
semilogy(EbNodB,BER_64QAM, 'DisplayName', '64-QAM')
grid on
legend('show')
legend('Location', 'southwest')
ylim([1e-20, 1])
yline(1e-6, '-.', 'DisplayName', '1e-6')
xlabel('Eb/N0 [dB]')
ylabel('BER')

berawgn(10.6, 'psk', 2, 'nondiff')
berawgn(10.6, 'psk', 4, 'nondiff')
berawgn(14, 'psk', 8, 'nondiff')
berawgn(14.4, 'qam', 16)
berawgn(18.7, 'qam', 64)
