#!/bin/bash

PYTHONPATH=src EGRIN2API_SETTINGS=test_settings.cfg coverage run test/api_test.py xml && coverage xml --include=src/app.py
