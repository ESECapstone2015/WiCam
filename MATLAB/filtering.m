
% 
% 
% % Generate a windowed sinc
% td_filter_len = 300;   
% td_filter_speed = 5;
% S = sinc(linspace(-td_filter_speed*pi,td_filter_speed*pi,td_filter_len+1));
% S = S / (sum(S));
% W = f.HammingWindow(0,td_filter_len,0);   
% S = S .* W;
% 
% % Plot windowed sinc
% S_t = f.to_time_domain(S, source_sample_rate, td_resolution, [-1,1], [-1,1]);
% gfx.figDiscreteAndTD('Windowed Sinc', S, S_t, [-inf,inf], [-inf,inf], source_plot_color);

% 
% S_t = f.to_time_domain(S, iq_sample_rate, td_resolution, source_sample_span, source_sample_span);
% 
% 
% figure(2);
% title = 'Sinc Window'; 
% set(gcf,'numbertitle','off','name',title);
% subplot(2, 1, 1);
% stem(S, 'Color',[0 0 0], 'MarkerFaceColor',[0 0 0]);
% subplot(2, 1, 2);
% plot(S);


%Y_t = conv(X_t(:,2), S);
%X_t(:,2) = Y_t(1:size(X_t,1));



% fc = 22.5e6;
% tau = 1/(2*pi*fc);
% a = 1/(iq_sample_rate*td_pts_per_sample*tau);

%X_t = filter(a, [1 a-1], X_t);

%figure(1);
%subplot(2, 1, 2);
%plot(X_t(:,1),X_t(:,2), 'Color',[0 0 0], 'LineWidth',1);
%axis([0,+inf,0,255]);