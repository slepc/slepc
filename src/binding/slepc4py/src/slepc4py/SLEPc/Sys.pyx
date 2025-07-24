# -----------------------------------------------------------------------------

cdef class Sys:
    """Sys."""

    @classmethod
    def getVersion(
        cls,
        devel: bool = False,
        date: bool = False,
        author: bool = False,
    ) -> tuple[int, int, int]:
        """Return SLEPc version information.

        Not collective.

        Parameters
        ----------
        devel
            Additionally, return whether using an in-development version.
        date
            Additionally, return date information.
        author
            Additionally, return author information.

        Returns
        -------
        major: int
            Major version number.
        minor: int
            Minor version number.
        micro: int
            Micro (or patch) version number.

        See Also
        --------
        slepc.SlepcGetVersion, slepc.SlepcGetVersionNumber

        """
        cdef char cversion[256]
        cdef PetscInt major=0, minor=0, micro=0, release=0
        CHKERR( SlepcGetVersion(cversion, sizeof(cversion)) )
        CHKERR( SlepcGetVersionNumber(&major, &minor, &micro, &release) )
        out = version = (toInt(major), toInt(minor), toInt(micro))
        if devel or date or author:
            out = [version]
            if devel:
                out.append(not <bint>release)
            if date:
                vstr = bytes2str(cversion)
                if release != 0:
                    date = vstr.split(",", 1)[-1].strip()
                else:
                    date = vstr.split("Git Date:")[-1].strip()
                out.append(date)
            if author:
                author = bytes2str(SLEPC_AUTHOR_INFO).split('\n')
                author = tuple([s.strip() for s in author if s])
                out.append(author)
        return tuple(out)

    @classmethod
    def getVersionInfo(cls) -> dict[str, bool | int | str]:
        """Return SLEPc version information.

        Not collective.

        Returns
        -------
        info: dict
            Dictionary with version information.

        See Also
        --------
        slepc.SlepcGetVersion, slepc.SlepcGetVersionNumber

        """
        version, dev, date, author = cls.getVersion(True, True, True)
        return dict(major      = version[0],
                    minor      = version[1],
                    subminor   = version[2],
                    release    = not dev,
                    date       = date,
                    authorinfo = author)

    # --- xxx ---

    @classmethod
    def isInitialized(cls) -> bool:
        """Return whether SLEPc has been initialized.

        Not collective.

        See Also
        --------
        isFinalized

        """
        return toBool(SlepcInitializeCalled)

    @classmethod
    def isFinalized(cls) -> bool:
        """Return whether SLEPc has been finalized.

        Not collective.

        See Also
        --------
        isInitialized

        """
        return toBool(SlepcFinalizeCalled)

    # --- xxx ---

    @classmethod
    def hasExternalPackage(cls, package: str) -> bool:
        """Return whether SLEPc has support for external package.

        Not collective.

        Parameters
        ----------
        package
            The external package name.

        See Also
        --------
        slepc.SlepcHasExternalPackage

        """
        cdef const char *cpackage = NULL
        package = str2bytes(package, &cpackage)
        cdef PetscBool has = PETSC_FALSE
        CHKERR( SlepcHasExternalPackage(cpackage, &has) )
        return toBool(has)

# -----------------------------------------------------------------------------
