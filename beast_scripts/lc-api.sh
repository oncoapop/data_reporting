#!/bin/bash

# Script to call Lab Collector APIs


GET
http://labcollector.online/apariciolab/webservice/v2/index.php?v=2&module=strains

X-LC-APP-Auth=ea85fe02636a422ff1e1850c8d78f88cc38c6f7cbd9b2c5959e42e2d6a76f6f4
Accept=text/xml 

#This gives a response



# API for slack as an example
#webhook_url="https://hooks.slack.com/services/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"

# This posts to my website channel on slack for realtime monitoring and for historical record
# Usage: slackpost "<webhook_url>" "<channel>" "<message>"


escapedText=$(echo $text | sed 's/"/\"/g' | sed "s/'/\'/g" | tr -d '"')
echo "Escaped text="$escapedText

json="{\"channel\": \"$channel\", \"text\": \"$escapedText\"}"
echo $json

curl -s -d "payload=$json" "$webhook_url"
