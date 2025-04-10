# Use the official Apache HTTPD image
FROM httpd:2.4

# Install Perl and mod_cgi to enable Perl CGI, and other stuff
RUN apt-get update && apt-get install -y \
    perl \
    libdbi-perl \
    libdbd-mysql-perl \
    libxml-simple-perl \
    libbio-perl-perl \
    libio-stringy-perl \
    libmath-vecstat-perl \
    build-essential \
    cpanminus \
    libapache2-mod-perl2 \
    cron \
    && rm -rf /var/lib/apt/lists/*

# Enable mod_cgi to process Perl scripts
#RUN a2enmod cgi
#RUN cpanm CGI
RUN cpanm CGI Email::Simple Email::Sender::Simple Email::Sender::Transport::SMTP

RUN cpanm Bio::Perl

RUN apt-get update && \
    apt-get install -y \
    procps && \
    apt-get clean

# Copy Perl CGI scripts into the container
COPY ./cgi-bin /usr/local/apache2/cgi-bin

# Set proper permissions for CGI scripts
RUN chmod +x /usr/local/apache2/cgi-bin/*

# Copy htdocs files into the container
COPY ./htdocs /usr/local/apache2/htdocs

# copy crontab file into the container
COPY crontab /etc/cron.d/crontab

# set permissions for the crontab file
RUN chmod 0644 /etc/cron.d/crontab

# Apply cron job
RUN crontab /etc/cron.d/crontab

# Set permissions for the CGI scripts to be executable
#RUN chmod -R 755 /usr/local/apache2/cgi-bin

# Set permissions for outputs and stat 
RUN chmod -R 777 /usr/local/apache2/htdocs/outputs
RUN chmod 777 /usr/local/apache2/htdocs/stats/webpssm.log

# Enable CGI module in Apache
RUN echo "LoadModule cgi_module modules/mod_cgi.so" >> /usr/local/apache2/conf/httpd.conf && \
    echo "AddHandler cgi-script .cgi .pl" >> /usr/local/apache2/conf/httpd.conf && \
    sed -i 's/DirectoryIndex index.html/DirectoryIndex index.html index.cgi/' /usr/local/apache2/conf/httpd.conf && \
    sed -i '/<Directory "\/usr\/local\/apache2\/htdocs">/,/<\/Directory>/ s/AllowOverride None/AllowOverride All/' /usr/local/apache2/conf/httpd.conf && \
    sed -i '/<Directory "\/usr\/local\/apache2\/cgi-bin">/,/<\/Directory>/ s/Options None/Options +ExecCGI/' /usr/local/apache2/conf/httpd.conf

# Expose port 80 for HTTP traffic
EXPOSE 80

# Start Apache in the foreground and enable cron service
CMD service cron start && httpd -D FOREGROUND