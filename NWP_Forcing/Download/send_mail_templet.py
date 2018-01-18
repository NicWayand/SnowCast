import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText

# TO USE: Rename this file to 'send_mail.py' and update your personal settings below

def send(body): 
    server = smtplib.SMTP('smtp.gmail.com', 587) # If using gmail
    server.starttls()
    server.login("senderEmail@gmail.com", "password")
    
    msg = MIMEMultipart()
    msg['From'] = "senderEmail@gmail.com.com"
    msg['To'] = "targetEmail@gmail.com"
    msg['Subject'] = "SnowCast Alert"
    
    msg.attach(MIMEText(body, 'plain'))
 
    server.sendmail("senderEmail@gmail.com", "targetEmail@gmail.com", str(msg))
    server.quit()
