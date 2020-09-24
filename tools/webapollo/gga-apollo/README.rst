Galaxy-apollo
=============

Local branch of `Galaxy Genome Annotation's apollo galaxy tools <https://github.com/galaxy-genome-annotation/galaxy-tools>`__

Galaxy tools to interface with Apollo.
Uses `python-apollo <https://github.com/galaxy-genome-annotation/python-apollo>`__ for most of it.

Environment
-----------

The following environment variables can be set:

+--------------------------------+-----------------------------------------------------------+
| ENV                            | Use                                                       |
+================================+===========================================================+
| ``$GALAXY_WEBAPOLLO_URL``      | The URL at which Apollo is accessible, internal to Galaxy |
|                                | and where the tools run. Must be absolute, with FQDN and  |
|                                | protocol.                                                 |
+--------------------------------+-----------------------------------------------------------+
| ``$GALAXY_WEBAPOLLO_USER``     | The admin user which Galaxy should use to talk to Apollo. |
|                                |                                                           |
+--------------------------------+-----------------------------------------------------------+
| ``$GALAXY_WEBAPOLLO_PASSWORD`` | The password for the admin user.                          |
+--------------------------------+-----------------------------------------------------------+
| ``$ARROW_GLOBAL_CONFIG_PATH``  | Path to a python-apollo/arrow conf file. Use in place of  |
|                                | ``$GALAXY_WEBAPOLLO_URL``, ``$GALAXY_WEBAPOLLO_USER``     |
|                                | and ``$GALAXY_WEBAPOLLO_PASSWORD``.                       |
+--------------------------------+-----------------------------------------------------------+
| ``$GALAXY_WEBAPOLLO_EXT_URL``  | May be relative or absolute.                              |
|                                | The external URL at which Apollo is accessible to end     |
|                                | users.                                                    |
+--------------------------------+-----------------------------------------------------------+
| ``$GALAXY_SHARED_DIR``         | Directory shared between Galaxy and Apollo, used to       |
|                                | exchange JBrowse instances. If not set, JBrowse data will |
|                                | be zipped and sent to the remote server.                  |
+--------------------------------+-----------------------------------------------------------+
| ``$GALAXY_APOLLO_ORG_SUFFIX``  | Set to 'id' if you want organism names to be suffixed     |
|                                | with user id to avoid name collisions. Set to 'email' to  |
|                                | use user email as suffix. Leave empty for no suffix.      |
+--------------------------------+-----------------------------------------------------------+

License
-------

All python scripts, wrappers, and the webapollo.py are licensed under
MIT license.
