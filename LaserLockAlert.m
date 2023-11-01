function LaserLockAlert
% User input
source = 'laserlocking@gmail.com';              %from address (gmail)
destination = 'yosri.benaicha@anu.edu.com';     %to address (any mail service)
% myEmailPassword = 'xgwo nnuu lvhs niok';        %the password to the 'from' account
subj = 'This is the subject line of the email'; % subject line
msg = 'This is the main body of the email.';    % main body of email.

% set up SMTP service for Gmail
setpref('Internet','E_mail',source);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username','laserlocking@gmail.com');
setpref('Internet','SMTP_Password','xgwo nnuu lvhs niok');

% Gmail server.
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

sendmail(destination,subj,msg);
% sendmail('Ryan.husband@anu.edu.au','laser locking status','Hurry up come back!!! LASER is unlocked');
% sendmail('Lorcan.conlon@anu.edu.au','laser locking status','Hurry up come back!!! LASER is unlocked');
% % sendmail('Zain.mehdi@anu.edu.au','laser locking status','Hurry up come back!!! LASER is unlocked');
% sendmail('Rhys.eagle@anu.edu.au','laser locking status','Hurry up come back!!! LASER is unlocked');
% sendmail('yosri.benaicha@anu.edu.au','laser locking status','Hurry up and come back!!! LASER is unlocked');
% sendmail('Samuel.legge@anu.edu.au','laser locking status','Hurry up come back!!! LASER is unlocked')

% Slack message
slackToken = 'xoxb-17496591987-5637707962160-Y0n0GlwNNU25pvFFBo03LUIp'; % filled with the token
channelID = 'C05J4DN718U'; % this to be filled with the channel ID
message = 'Hurry up and come back!!! Gravy Laser is unlocked';

% Create a structure for the payload
payload = struct('channel', channelID, 'text', message);

% Create HTTP request
header(1) = matlab.net.http.field.ContentTypeField('application/json;charset=UTF-8');
header(2) = matlab.net.http.field.AuthorizationField('Authorization', ['Bearer ', slackToken]);
request = matlab.net.http.RequestMessage('POST', header, payload); % let MATLAB handle the encoding

% Send HTTP request
uri = matlab.net.URI('https://slack.com/api/chat.postMessage');
response = request.send(uri);

% Cheking the respons body is ok! (1)
% if response.StatusCode == matlab.net.http.StatusCode.OK
%     disp('Message sent successfully.')
%     disp('Response body:')
%     disp(response.Body.Data)
% else
%     disp('Failed to send message:')
%     disp(response.Body)
% end
 

end
