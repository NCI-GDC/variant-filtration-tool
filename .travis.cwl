language: python

python:
  - '3.5'

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
  - secure: Uo8xIwl+9SNy8UzIFn8J7+XpPWCX6jcnCYu2cAZq51QWEH75kBfrlE89u2Hu3i5eFbE/ufhTjXme5FcJpBfaW1Msa8lFt0u8UHdlvPcpP7Rbo5D+EqVJd3zviYa92s8qK4LWhjo6ZWcCISBa8zG8Ki3lCEob8b/mlnlYfL2PI2pqLmZA+JMa3LoGV5U/gBM1i7yvSu8xXq28aiugRpovnn0cT9QH7pTXii+HfTbCt4gSsNJaU8IDMWY3EK2vlmE9CT8TbO/TQNxgxbPHTjq7JInau9xJT+rD98pHVsxoW6gAV2J41GzfOKN/s0XSbGZhKVohNHf36akigdRPw0KXmlKcxQKVkDOcvKUA7Z/C1JmWcTbARa/Hmpwlw4CzP16O9aTpVbSGE0lYJjUkeQRceScKA9R1nCzk9XbzdtiQ4NDJzOl8LSungh6Jbwh2aWeJa7uBve0m9Zb37yZhM2UFY96WIvLG/uvUwQNHvJTbELbrjrixoym69nHuzJzL1AiYkvf7T9+MWS8118Fve37liZN8yZA8AhR7aL3Fx96wqM07dovM/xUg01evzrCcVIlQwTE2GqDjpqmC2Z67jlo4hv9PPOqQGO0dVuAZo7Wb1VN9VhHj3gEGEQS8YOOV2JCWtVDfUm2QSliO0T3ViKIbuqBIfQNdp5XjBwKpHMzNNpU=
  - secure: bD9LcjAhAN1emv8nJrIaBTSxOJ94N1kP9ZUUgqcGg6QmKPHKmjS1uZo5RvkNPrG843FHE5qT2cpViFzNip28XIqr3RAsMlTzxJWNdPxt/WBbygDLmKEvOPd+cGad7jrIZKKV5IVEXRx7tvfp9QNP01cwhh0/HovqPVtyL/bBSaMjNvlgYlWRaw36ymNDifmhUFQmHYllqXMveWDoaIkyyuUafUITATD9BJQ/KMV86wxNv1JfiRkYkVb4wT4ZS9RyOtuvmfbIw5pPuc/5KtvPAJghBmnMoCcuYN13t29cUMqhUlA4+RLruItUQTCJU0DiLiBBIjde5W9p/O/HrpK2z6dTBS7OIjsMCRG5V348dW6D/f+cUXHrE8mPWLoq4UCgymZPSGEgNBj7qxujhM8pnfizlIoW5/DKAfmmSERd4JKTacTv8vZ8kE+xcKedTLpfNcgWpuncU8rbSqtc7lyv+C6ZGAW2m/MrcedI8Ao/B7ju0pR5q8DQWdziboYamz62Ix5mQZLUtMc6HLEaYncJ3RFtAT9T6VupEaFpM8zyygDXPG3bznYiAQg3kZw+7F7P/DYx3tjlDq9tkfbJUpHgQbWuHyaI+6FwqgNtKYEnHg/W3WBwylCJmjxIlwKGcdIk4ohIfLsOx10EE/4nxs+bG5stQyJQCNisQ9W86GO0e6k=
