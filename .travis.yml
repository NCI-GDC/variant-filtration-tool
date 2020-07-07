language: python

python:
  - '3.6'

services: docker

git:
  depth: false

install: make build

stages:
- name: build
- name: publish
  if: branch = master
- name: publish-staging
  if: branch != master
jobs:
  include:
  - stage: build
    name: Build
    script:
    - make version version-pypi
    - make test-docker
  - stage: publish-staging
    name: Publish staging image
    script:
    - docker login -u="$QUAY_USERNAME" -p="$QUAY_PASSWORD" quay.io
    - make publish-staging
  - stage: publish
    name: Publish production image
    script:
    - docker login -u="$QUAY_USERNAME" -p="$QUAY_PASSWORD" quay.io
    - make publish-release
env:
  global:
  - BRANCH_NAME=${TRAVIS_PULL_REQUEST_BRANCH:-$TRAVIS_BRANCH}
