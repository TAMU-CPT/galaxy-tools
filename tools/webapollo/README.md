# galaxy-apollo

Galaxy tools to interface with Apollo
The webapollo.py file is also [separately available](https://github.com/erasche/python-apollo) as a pip-installable package.

## ENV

The following environment variables must be set:


ENV                           | Use
---                           | ---
 `$GALAXY_WEBAPOLLO_URL`      | The URL at which Apollo is accessible, internal to Galaxy and where the tools run. Must be absolute, with FQDN and protocol.
 `$GALAXY_WEBAPOLLO_USER`     | The admin user which Galaxy should use to talk to Apollo.
 `$GALAXY_WEBAPOLLO_PASSWORD` | The password for the admin user.
 `$GALAXY_WEBAPOLLO_EXT_URL`  | The external URL at which Apollo is accessible to end users. May be relative or absolute.
 `$GALAXY_SHARED_DIR`         | Directory shared between Galaxy and Apollo, used to exchange JBrowse instances.

## License

All python scripts, wrappers, and the webapollo.py are licensed under MIT license.
