from django import forms

class ModelForm(forms.Form):
    type_select = forms.ChoiceField(choices = ( (0, 'Single Strategy [Click for List]'),
                                                (1, 'All Strategies') ), widget=forms.RadioSelect() )
    int_select = forms.ChoiceField(choices=( (1, '1. Xpert for smear-positive'),
                                (2, '2. Xpert for HIV+'),
                                (3, '3. Xpert for previously treated'),
                                (4, '4. Xpert for sm-neg HIV+ or prev tx'),
                                (5, '5. Xpert for all HIV+ or prev tx'),
                                (6, '6. Xpert for smear-negative'),
                                (7, '7. Xpert for all'),
                                (8, '8. Xpert for all, same-day') ), widget=forms.RadioSelect() )
    t_inc = forms.FloatField(
            label='Target TB incidence, per 100,000',
            max_value=100000, min_value=0.0, initial=250, widget=forms.TextInput(attrs={'size': '5'}))
    t_mdr = forms.FloatField(
            label='Target MDR-TB prevalence among new cases, %',
            max_value=100.0, min_value=0.0, initial=3.7, widget=forms.TextInput(attrs={'size': '5'}))
    t_hiv = forms.FloatField(
            label='Target adult HIV prevalence, %',
            max_value=100.0, min_value=0.0, initial=0.83, widget=forms.TextInput(attrs={'size': '5'}))
    t_drug1_cost = forms.FloatField(
            label='Treatment of one patient with first-line drugs, $',
            min_value=0.0, initial=500, widget=forms.TextInput(attrs={'size': '5'}))
    t_drug2_cost = forms.FloatField(
            label ='Treatment of one patient with retreatment ("category 2") regimen, $', min_value=0.0,
            initial=1000, widget=forms.TextInput(attrs={'size':'5'}))
    t_drug3_cost = forms.FloatField(
            label ='Treatment of one patient with second-line (MDR) drugs, $', min_value=0.0,
            initial=5000, widget=forms.TextInput(attrs={'size':'5'}))
    t_outpt_cost = forms.FloatField(
            label ='One outpatient visit (e.g., for TB diagnosis), $', min_value=0.0,
            initial=10, widget=forms.TextInput(attrs={'size':'2'}))
    t_sm_cost = forms.FloatField(
            label='Full sputum smear evaluation (e..g, collection & evaluation of 2 smears), $', min_value=0.0, initial=2,
            widget=forms.TextInput(attrs={'size':'2'}))
    t_gxp_cost = forms.FloatField(
            label='Single Xpert MTB/RIF test, $', min_value=0.0, initial=15,
            widget=forms.TextInput(attrs={'size':'2'}))
    t_sdgxp_cost = forms.FloatField(
            label='Single Xpert, including extra costs to make results available same day, $', min_value=0.0,
            initial=30, widget=forms.TextInput(attrs={'size':'2'}))


