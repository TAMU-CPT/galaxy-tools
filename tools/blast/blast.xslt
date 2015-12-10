<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

    <xsl:output method="xml" encoding="UTF-8" indent="yes"/>

    <xsl:template match="@*|node()">
        <xsl:copy>
            <xsl:apply-templates select="@*|node()"/>
        </xsl:copy>
    </xsl:template>

    <xsl:template match="Iteration_hits/Hit">
        <xsl:if test='Hit_hsps[Hsp[number(Hsp_evalue) &lt; 1e-07 and number(Hsp_num) &lt; 3]]'>
            <xsl:copy>
                <xsl:apply-templates select="@*|node()"/>
            </xsl:copy>
        </xsl:if>
    </xsl:template>

    <xsl:template match="Iteration_hits/Hit/Hit_hsps/Hsp">
        <xsl:if test='self::*[number(Hsp_evalue) &lt; 1e-07 and number(Hsp_num) &lt; 3]'>
            <xsl:copy>
                <xsl:apply-templates select="@*|node()"/>
            </xsl:copy>
        </xsl:if>
    </xsl:template>

</xsl:stylesheet>
