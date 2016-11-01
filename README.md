# cmonkey2 API (cmapi)

## Description

cmapi is a project to provide the data of a cmonkey2 run or EGRIN2 run through
an abstract web service API.

## Technical concept

The API is implemented as a Flask application with SQL Alchemy implementing
a database independent layer on top of any supported relational database.
A converter script helps importing the sqlite database into a database
of the user's choice.

## Frontend integration

The API offers various formats as exports: JSON, CSV, HTML, SVG and PNG
depending on the request.

A Javascript API allows the configuration of the returned components
for rendering in the browser.

There will also be a Wordpress plugin to support users in configuring
the components through the Wordpress administration interface

