from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
import smtplib


def notify(subject, body, attachment):
    """Send a notification email to the CPT's Galaxy admin.

    Can be useful for finding odd edge cases. TODO: replace with proper logging/alerting.

    subject and body are just strings. Attachment should be a dictionary with a
    'content' and a 'filename'
    """
    msg = MIMEMultipart()
    msg['Subject'] = subject
    msg['From'] = 'galaxy@cpt.tamu.edu'
    msg['To'] = 'esr@tamu.edu'

    msg.attach(MIMEText(body))

    attachment = MIMEText(attachment['content'])
    attachment.add_header('Content-Disposition', 'attachment', filename=attachment['filename'])
    msg.attach(attachment)

    # to send
    mailer = smtplib.SMTP()
    mailer.connect()
    mailer.sendmail("root@cpt.tamu.edu", "esr@tamu.edu", msg.as_string())
    mailer.close()
