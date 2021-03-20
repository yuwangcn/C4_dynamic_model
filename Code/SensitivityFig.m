function SensitivityFig(Timex,RCC)
figure;%Fig 5
subplot(1,3,1);
plot(Timex,RCC(:,1),'k.');
title('ki\_gs');xlim([0,1800]);ylim([0,1.5]);
xlabel('Time (s)');
ylabel('SC' );
subplot(1,3,2);
plot(Timex,RCC(:,2),'k.');
title('1/TaoRubisco');xlim([0,1800]);ylim([0,1.5]);
xlabel('Time (s)');
ylabel('SC' );
subplot(1,3,3);
plot(Timex,RCC(:,3),'k.');
title('RP');xlim([0,1800]);ylim([0,1.5]);
xlabel('Time (s)');
ylabel('SC' );

figure;%Fig 6-1
subplot(2,3,1);
plot(Timex,RCC(:,7),'k.');
title('PEPC');xlim([0,1800]);ylim([0,1.5]);
xlabel('Time (s)');
ylabel('SC' );
subplot(2,3,2);
plot(Timex,RCC(:,9),'k.');
title('PPDK');xlim([0,1800]);ylim([0,1.5]);
xlabel('Time (s)');
ylabel('SC' );
subplot(2,3,3);
plot(Timex,RCC(:,10),'k.');
title('MDH');xlim([0,1800]);ylim([0,1.5]);
xlabel('Time (s)');
ylabel('SC' );
subplot(2,3,4);
plot(Timex,RCC(:,11),'k.');
title('ME');xlim([0,1800]);ylim([0,1.5]);
xlabel('Time (s)');
ylabel('SC' );
subplot(2,3,5);
plot(Timex,RCC(:,17),'k.');
title('Mutase&Enolase');xlim([0,1800]);ylim([0,1.5]);
xlabel('Time (s)');
ylabel('SC' );

figure;%Fig 6-2
subplot(2,3,1);
plot(Timex,RCC(:,4),'k.');
title('Rubisco');xlim([0,1800]);ylim([0,1.5]);
xlabel('Time (s)');
ylabel('SC' );
subplot(2,3,2);
plot(Timex,RCC(:,12),'k.');
title('DAPDH');xlim([0,1800]);ylim([0,1.5]);
xlabel('Time (s)');
ylabel('SC' );
subplot(2,3,3);
plot(Timex,RCC(:,13),'k.');
title('SBPase');xlim([0,1800]);ylim([0,1.5]);
xlabel('Time (s)');
ylabel('SC' );
subplot(2,3,4);
plot(Timex,RCC(:,14),'k.');
title('FBPase');xlim([0,1800]);ylim([0,1.5]);
xlabel('Time (s)');
ylabel('SC' );
subplot(2,3,5);
plot(Timex,RCC(:,15),'k.');
title('PRK');xlim([0,1800]);ylim([0,1.5]);
xlabel('Time (s)');
ylabel('SC' );

